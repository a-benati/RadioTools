#!/usr/bin/env python3

import argparse
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np

def crop_fits_using_reference(big_fits, small_fits):
    with fits.open(big_fits) as big_hdul, fits.open(small_fits) as small_hdul:
        big_hdu, small_hdu = big_hdul[0], small_hdul[0]
        
        big_wcs = WCS(big_hdu.header).celestial
        small_wcs = WCS(small_hdu.header).celestial
        
        big_data = big_hdu.data.squeeze()
        small_data = small_hdu.data.squeeze()
        
        big_pixel_scale = np.sqrt(np.abs(big_wcs.pixel_scale_matrix[0, 0] * big_wcs.pixel_scale_matrix[1, 1]))
        small_pixel_scale = np.sqrt(np.abs(small_wcs.pixel_scale_matrix[0, 0] * small_wcs.pixel_scale_matrix[1, 1]))
        print(f"Big pixel scale: {big_pixel_scale}")
        print(f"Small pixel scale: {small_pixel_scale}")
        
        scale_factor = small_pixel_scale / big_pixel_scale
        
        ny, nx = small_data.shape  # Dimensions of small image
        print(f"Original small image dimensions: {nx} x {ny}")
        
        scaled_nx, scaled_ny = int(nx / scale_factor), int(ny / scale_factor)
        print(f"Scaled small image dimensions: {scaled_nx} x {scaled_ny}")
        
        # Get the corners of the small image in pixel coordinates of the big image
        corners_small = np.array([[0, 0], [nx-1, 0], [0, ny-1], [nx-1, ny-1]])
        corners_sky = small_wcs.pixel_to_world(corners_small[:, 0], corners_small[:, 1])
        corners_big = big_wcs.world_to_pixel(corners_sky)
        
        xmin, xmax = int(np.min(corners_big[0])), int(np.max(corners_big[0]))
        ymin, ymax = int(np.min(corners_big[1])), int(np.max(corners_big[1]))
        
        # Ensure boundaries are within the image dimensions
        xmin, xmax = max(0, xmin), min(big_data.shape[1], xmax)
        ymin, ymax = max(0, ymin), min(big_data.shape[0], ymax)
        print(f"Bounding box: ({xmin}, {ymin}) to ({xmax}, {ymax})")
        
        cropped_data = big_data[ymin+2:ymax, xmin+2:xmax]
        print(f"Cropped image dimensions: {cropped_data.shape}")
        
        new_header = big_hdu.header.copy()
        new_header['CRPIX1'] -= xmin
        new_header['CRPIX2'] -= ymin
        
        output_file = big_fits.replace(".fits", "-crop.fits")
        fits.writeto(output_file, cropped_data, new_header, overwrite=True)
        print(f"Cropped image saved as: {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Crop a FITS image using another smaller FITS as reference.")
    parser.add_argument("big_fits", help="Big fits image (to crop)")
    parser.add_argument("small_fits", help="Small fits image (reference for cropping)")
    args = parser.parse_args()
    
    crop_fits_using_reference(args.big_fits, args.small_fits)
