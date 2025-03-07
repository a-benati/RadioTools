#!/usr/bin/env python3

import numpy as np
from astropy.io import fits
import subprocess
import argparse
import os
from casacore.tables import table

def flux_density_within_region(fits_file, region_files, reg2fits_path):
    # Open the FITS file and extract the image data and header
    hdu = fits.open(fits_file)[0]
    data = hdu.data.squeeze()  # Assuming 2D data
    header = hdu.header

    # Get pixel scale (in degrees)
    cdelt1 = np.abs(header['CDELT1'])  # Pixel scale in x (deg/pixel)
    cdelt2 = np.abs(header['CDELT2'])  # Pixel scale in y (deg/pixel)

    print(f"Pixel scale (CDELT1, CDELT2): {cdelt1}, {cdelt2}")

    # Compute pixel area in steradians
    pixel_area_sr = cdelt1 * cdelt2 * (np.pi / 180.0) ** 2
    print(f"Pixel area in steradians: {pixel_area_sr}")

    # Get beam major and minor axes (in degrees)
    bmaj = header['BMAJ']  # Major axis (degrees)
    bmin = header['BMIN']  # Minor axis (degrees)

    print(f"Beam major (BMAJ) and minor (BMIN): {bmaj}, {bmin}")

    # Compute beam area in steradians
    beam_area_sr = np.pi * (bmaj / 2) * (bmin / 2) * (np.pi / 180.0) ** 2
    print(f"Beam area in steradians: {beam_area_sr}")

    # Loop through each region file
    for region_file in region_files:
        # Call the external reg2fits.py script to create a mask FITS file
        mask_fits = os.path.splitext(fits_file)[0] + '_mask.fits'
        subprocess.run(['python3', reg2fits_path, region_file, fits_file])

        # Open the mask FITS file and extract the mask data
        mask_hdu = fits.open(mask_fits)[0]
        mask = mask_hdu.data.squeeze()  # Assuming 2D mask data

        # Apply the mask to the data to select the region
        selected_region_data = data[mask == 1]

        # Compute sum, mean, std, min, and max in Jy/beam
        pixel_sum_jy_per_beam = np.sum(selected_region_data)
        pixel_mean_jy_per_beam = np.mean(selected_region_data)
        pixel_std_jy_per_beam = np.std(selected_region_data)
        pixel_min_jy_per_beam = np.min(selected_region_data)
        pixel_max_jy_per_beam = np.max(selected_region_data)

        # Print statistics for the selected region
        print(f"\nRegion: {region_file}")
        print(f"Sum of pixels in region (Jy/beam): {pixel_sum_jy_per_beam}")
        print(f"Mean pixel value in region (Jy/beam): {pixel_mean_jy_per_beam}")
        print(f"Standard deviation of pixel values (Jy/beam): {pixel_std_jy_per_beam}")
        print(f"Minimum pixel value in region (Jy/beam): {pixel_min_jy_per_beam}")
        print(f"Maximum pixel value in region (Jy/beam): {pixel_max_jy_per_beam}")

        # Convert the sum from Jy/beam to Jy
        flux_density_jy = pixel_sum_jy_per_beam * (pixel_area_sr / beam_area_sr)
        print(f"Flux density in region (Jy): {flux_density_jy}")

        # Calculate the flux density with Casacore
        #ms_file = os.path.splitext(fits_file)[0] + '.ms'  # Assuming you have an MS file
        # ms_file = '/data/abell_3667/msfiles/1685906777_sdp_l0-cal.ms'
        # tb = table(ms_file)
        
        # Get the flux density using Casacore
        # You may need to adjust the following based on the structure of your MS
        # Example: Compute total flux density by summing the fluxes in the selected region
        # tb = table(ms_file, readonly=True)
        # data_column = tb.getcol('CORRECTED_DATA')  # Adjust this if necessary
        # if np.iscomplexobj(data_column):
        #     flux_sum = np.sum(data_column.real)  # Only take the real part
        # else:
        #     flux_sum = np.sum(data_column)

        # flux_density_casacore = flux_sum * (pixel_area_sr / beam_area_sr)
        
        # print(f"Flux density (using Casacore) in region (Jy): {flux_density_casacore}")

        # tb.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate flux density within a region.')
    parser.add_argument('-img', required=True, help='FITS image file')
    parser.add_argument('-reg', nargs='+', required=True, help='Region files (in .reg format)')
    parser.add_argument('--reg2fits_path', default='/packages/scripts/reg2fits.py',
                        help='Path to reg2fits.py script')

    args = parser.parse_args()
    flux_density_within_region(args.img, args.reg, args.reg2fits_path)
