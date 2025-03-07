#!/usr/bin/env python3

import argparse
import numpy as np
from astropy.io import fits
import os
from subprocess import call


def compute_polarisation_fraction(I, Q, U):
    """Compute the polarisation fraction as sqrt(Q^2 + U^2) / I"""
    P_lin = np.sqrt(Q**2 + U**2)
    frac_pol = P_lin / I if I != 0 else 0
    return frac_pol


def compute_polarisation_angle(Q, U):
    """Compute the polarisation angle in degrees as 0.5 * arctan(U / Q)"""
    pol_angle = 0.5 * np.arctan2(U, Q)  # Arctan2 handles the quadrant correctly
    pol_angle_deg = np.degrees(pol_angle)
    return pol_angle_deg


def load_fits_data(file_path):
    """Load FITS file data."""
    with fits.open(file_path) as hdulist:
        data = hdulist[0].data
    return data


def integrate_flux(image_data, mask_data):
    """Integrate flux within the region mask"""
    return np.sum(image_data * mask_data)


def apply_region_mask(region_file, image_file):
    """Convert a region file to a mask and apply it to an image"""
    mask_file = os.path.splitext(image_file)[0] + '_mask.fits'
    call(["python3", "/RadioTools/reg2fits.py", "-reg", region_file, "-img", image_file, "-o", mask_file])
    
    mask_data = load_fits_data(mask_file)
    return mask_data


def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Compute integrated polarisation fraction and angle for multiple channels")
    parser.add_argument("--prefix", required=True, help="Prefix for the channel FITS files (e.g., '/data/abell_3667/polcal/3c286/3c286_CORRECTED')")
    parser.add_argument("--region", required=True, help="Region file to select the source area")
    parser.add_argument("--num_channels", type=int, default=15, help="Number of channels to process (default: 15)")
    parser.add_argument("--output", required=True, help="Output text file to store the results")

    args = parser.parse_args()

    # Open the output file for writing
    with open(args.output, 'w') as f:
        f.write("Channel\tPolarisation Fraction\tPolarisation Angle (deg)\n")
        
        for channel in range(args.num_channels):
            # Create the filenames for I, Q, U, V
            I_file = f"{args.prefix}-{'{:04d}'.format(channel)}-I-image_pbcor.fits"
            Q_file = f"{args.prefix}-{'{:04d}'.format(channel)}-Q-image_pbcor.fits"
            U_file = f"{args.prefix}-{'{:04d}'.format(channel)}-U-image_pbcor.fits"
            V_file = f"{args.prefix}-{'{:04d}'.format(channel)}-V-image_pbcor.fits"  # Not used but listed for completeness
            
            # Load the Stokes parameters
            I_data = load_fits_data(I_file)
            Q_data = load_fits_data(Q_file)
            U_data = load_fits_data(U_file)

            # Apply the region mask
            mask_data = apply_region_mask(args.region, I_file)

            # Integrate the flux inside the region for I, Q, U
            I_flux = integrate_flux(I_data, mask_data)
            Q_flux = integrate_flux(Q_data, mask_data)
            U_flux = integrate_flux(U_data, mask_data)

            # Calculate the polarisation fraction and angle
            pol_fraction = compute_polarisation_fraction(I_flux, Q_flux, U_flux)
            pol_angle = compute_polarisation_angle(Q_flux, U_flux)

            # Write the results for this channel to the file
            f.write(f"{channel:04d}\t{pol_fraction:.6f}\t{pol_angle:.6f}\n")
            print(f"Processed channel {channel:04d}")

    print(f"Results saved to {args.output}")


if __name__ == "__main__":
    main()
