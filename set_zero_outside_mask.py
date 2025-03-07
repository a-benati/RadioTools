#!/usr/bin/env python3

import argparse
from astropy.io import fits
import numpy as np

# Define the function to apply the mask and zero the pixels outside of it
def apply_mask(image_file, mask_file, output_file):
    # Load the image and the mask
    with fits.open(image_file) as hdul:
        image_data = hdul[0].data

    with fits.open(mask_file) as hdul:
        mask_data = hdul[0].data

    # Apply the mask: zero the pixels outside the mask
    masked_image_data = np.where(mask_data > 0, image_data, 0)

    # Save the new masked image
    hdu = fits.PrimaryHDU(masked_image_data)
    hdu.writeto(output_file, overwrite=True)
    print(f"Output image saved as {output_file}")

# Main function to handle command-line arguments
def main():
    parser = argparse.ArgumentParser(description="Apply a mask to an image, zeroing pixels outside the mask")
    parser.add_argument("--image", help="Path to the FITS image file")
    parser.add_argument("--mask", help="Path to the FITS mask file")
    parser.add_argument("--output", help="Path to save the masked FITS image (default: 'zeros_outside_mask_image.fits')", default="zeros_outside_mask_image.fits")

    args = parser.parse_args()

    # Apply the mask to the image
    apply_mask(args.image, args.mask, args.output)

if __name__ == "__main__":
    main()
