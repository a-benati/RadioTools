#!/usr/bin/env python3

import argparse
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.wcs import WCS

def read_region_file(region_file, wcs):
    # Read the region file and create a mask for the ellipse
    with open(region_file, 'r') as file:
        lines = file.readlines()
    
    # Extract the parameters from the region file
    for line in lines:
        if line.startswith('ellipse'):
            # Example line: ellipse(-157.215417228, 30.509140655, 20.5128", 8.6509", 82.430084)
            parts = line.split('(')[1].split(')')[0].split(',')
            ra = float(parts[0].strip())
            dec = float(parts[1].strip())
            semi_major = float(parts[2].strip().replace('"', '').strip()) * u.arcsec
            semi_minor = float(parts[3].strip().replace('"', '').strip()) * u.arcsec
            angle = float(parts[4].strip())

            # Convert the RA/Dec to pixel coordinates
            center = SkyCoord(ra, dec, unit='deg', frame='fk5')
            
            # Ensure we are passing the correct number of dimensions
            world_coords = np.array([[center.ra.deg, center.dec.deg]])  # Wrap in a 2D array
            x_center, y_center = wcs.world_to_pixel(world_coords)  # Call with 2D array
            
            x_center = x_center[0]
            y_center = y_center[0]

            # Create an array to hold the mask
            mask = np.zeros(wcs.array_shape, dtype=bool)

            # Create a grid of pixel coordinates
            y_indices, x_indices = np.indices(mask.shape)

            # Calculate the distance from the center
            a = semi_major.to(u.pix, wcs).value  # Semi-major axis in pixels
            b = semi_minor.to(u.pix, wcs).value   # Semi-minor axis in pixels

            # Ellipse equation
            ellipse_eq = ((x_indices - x_center) * np.cos(np.radians(angle)) +
                           (y_indices - y_center) * np.sin(np.radians(angle)))**2 / a**2 + \
                         ((x_indices - x_center) * np.sin(np.radians(angle)) -
                           (y_indices - y_center) * np.cos(np.radians(angle)))**2 / b**2

            # Set the mask to True for points inside the ellipse
            mask[ellipse_eq <= 1] = True
            return mask
    return None

def calculate_polarization(I_file, Q_file, U_file, region_file):
    # Load the Stokes I, Q, and U images
    I_data = fits.getdata(I_file)
    Q_data = fits.getdata(Q_file)
    U_data = fits.getdata(U_file)

    # Create WCS for the images
    header = fits.getheader(I_file)
    wcs = WCS(header)

    # Create a mask based on the region file
    mask = read_region_file(region_file, wcs)

    if mask is None:
        print("Region mask could not be created.")
        return

    # Integrated Stokes parameters within the masked region
    I_integrated = np.sum(I_data[mask])
    Q_integrated = np.sum(Q_data[mask])
    U_integrated = np.sum(U_data[mask])

    # Calculate the polarization fraction
    P = np.sqrt(Q_integrated**2 + U_integrated**2) / I_integrated if I_integrated != 0 else 0

    # Calculate the polarization angle
    theta = 0.5 * np.arctan2(U_integrated, Q_integrated) * (180 / np.pi)  # Convert to degrees

    print(f'Integrated Stokes I: {I_integrated:.2f}')
    print(f'Integrated Stokes Q: {Q_integrated:.2f}')
    print(f'Integrated Stokes U: {U_integrated:.2f}')
    print(f'Polarization Fraction: {P:.2f}')
    print(f'Polarization Angle: {theta:.2f} degrees')

    # Expected values
    expected_P = 0.10
    expected_theta = 33.0

    if np.isclose(P, expected_P, atol=0.01):
        print("Polarization fraction is as expected.")
    else:
        print("Polarization fraction is NOT as expected.")

    if np.isclose(theta, expected_theta, atol=1.0):
        print("Polarization angle is as expected.")
    else:
        print("Polarization angle is NOT as expected.")

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Calculate polarization fraction and angle from Stokes parameters.')
    parser.add_argument('-I', required=True, help='FITS file for Stokes I')
    parser.add_argument('-Q', required=True, help='FITS file for Stokes Q')
    parser.add_argument('-U', required=True, help='FITS file for Stokes U')
    parser.add_argument('-reg', required=True, help='Region file for masking')

    args = parser.parse_args()

    # Call the calculation function
    calculate_polarization(args.I, args.Q, args.U, args.reg)

if __name__ == '__main__':
    main()
