#!/usr/bin/env python3

import numpy as np
import argparse
import astropy.io.fits as fits
from scipy.ndimage import gaussian_filter
from skimage import measure

def create_contour_region(image_path, sigma_level, output_path):
    # Open the FITS file
    with fits.open(image_path) as hdul:
        # Assuming the data is in the primary HDU
        data = hdul[0].data

    # Calculate the contour level based on the mean and std deviation
    mean_val = np.mean(data)
    std_val = np.std(data)
    contour_level = mean_val + sigma_level * std_val

    # Smoothing the data (if necessary)
    smoothed_data = gaussian_filter(data, sigma=1)  # You can adjust the sigma here

    # Find contours at the specified level
    contours = measure.find_contours(smoothed_data[0], contour_level)  # Use [0] to take the first slice

    # Create a region file
    with open(output_path, 'w') as reg_file:
        reg_file.write("# Region file format: DS9 version 4.1\n")
        reg_file.write("global color=green\n")
        reg_file.write("image\n")

        for contour in contours:
            # Write each contour as a polygon
            contour = contour * 1.0  # Ensuring it's a float
            reg_file.write("polygon(")
            # Format the contour points for DS9
            reg_file.write(', '.join(f"{x},{y}" for x, y in contour))
            reg_file.write(")\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create a DS9 region file from contour data.')
    parser.add_argument('-i', '--image', required=True, help='Path to the input FITS file')
    parser.add_argument('-s', '--sigma', type=float, required=True, help='Sigma level for contour')
    parser.add_argument('-o', '--output', required=True, help='Path to the output region file')
    
    args = parser.parse_args()
    create_contour_region(args.image, args.sigma, args.output)
