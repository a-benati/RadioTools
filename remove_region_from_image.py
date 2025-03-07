#!/usr/bin/env python3

import argparse
import numpy as np
from astropy.io import fits
import os
import subprocess

def apply_mask(image_data, mask_file):
    """
    Apply the mask to the image data, setting pixels in the mask (1) to zero.
    
    Parameters:
    - image_data: numpy array of the image to be modified.
    - mask_file: path to the mask FITS file.
    
    Returns:
    - Modified image data with masked regions set to 0.
    """
    with fits.open(mask_file) as mask_hdu:
        mask_data = mask_hdu[0].data

        # Squeeze the image to remove any singleton dimensions (like (1, 10000, 10000))
        image_data_squeezed = np.squeeze(image_data)

        # Ensure the image and mask have the same shape
        if mask_data.shape != image_data_squeezed.shape:
            raise ValueError(f"Mask shape {mask_data.shape} does not match image shape {image_data_squeezed.shape}")

        # Set pixels where the mask is 1 to zero
        image_data_squeezed[mask_data == 1] = 0

        # Return the modified image, expanded back to original dimensions if necessary
        return image_data_squeezed

def convert_region_to_mask(region_file, image_file, mask_file):
    """
    Use the reg2fits.py script to convert a region file to a mask FITS file.
    
    Parameters:
    - region_file: Path to the region file.
    - image_file: Path to the FITS image file.
    - mask_file: Path to save the generated mask FITS file.
    """
    subprocess.run(['./reg2fits.py', '-reg', region_file, '-img', image_file, '-o', mask_file])

def main(image, output_image, region_files):
    # Open the input image (FITS)
    with fits.open(image) as hdulist:
        image_data = hdulist[0].data
        
        # Apply regions to the image data
        for region_file in region_files:
            # Generate a mask for each region
            mask_file = os.path.splitext(region_file)[0] + "_mask.fits"
            convert_region_to_mask(region_file, image, mask_file)
            
            # Apply the mask to the image data
            image_data = apply_mask(image_data, mask_file)
        
        # Determine the output file name
        if output_image:
            output_path = output_image
        else:
            base, ext = os.path.splitext(image)
            output_path = f"{base}_regions_removed.fits"
        
        # Save the modified image
        fits.writeto(output_path, image_data, hdulist[0].header, overwrite=True)
        print(f"Modified image saved to: {output_path}")

if __name__ == "__main__":
    # Set up argparse to handle input
    parser = argparse.ArgumentParser(description="Set pixels inside regions to 0 in the image")
    parser.add_argument("-r", "--region", nargs='+', required=True, help="Paths to region files to apply")
    parser.add_argument("-o", "--output_image", required=False, help="Path to save the modified output FITS image (optional)")
    parser.add_argument("-i", "--image", help="Path to the FITS image")

    args = parser.parse_args()

    # Run the main function
    main(args.image, args.output_image, args.region)
