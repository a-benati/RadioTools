#!/usr/bin/env python

from astropy.io import fits as pyfits
import argparse
import numpy as np

def combine_masks(output_mask, input_masks, invert_masks):
    # Determine the output mask file path
    if output_mask is None:
        output_mask = "mask.fits"
        
    # Initialize data_comb to None, will assign the first mask's data initially
    data_comb = None

    # Open each input mask file and combine using the "0 overwrites 1" logic
    for i, input_mask in enumerate(input_masks):
        with pyfits.open(input_mask) as fits:
            data = fits[0].data

            # Remove dimensions of size 1 (squeeze)
            data = np.squeeze(data)

            # Invert the mask if specified
            if invert_masks[i]:
                data = 1.0 - data  # Invert: 1 -> 0 and 0 -> 1

            if data_comb is None:
                # Start with the first mask data
                data_comb = np.ones_like(data)  # Start with all 1's for combination

            if data.shape != data_comb.shape:
                print(f"Warning: Shapes of input mask {input_mask} and output mask {output_mask} are different:")
                print(f"Shape of input mask: {data.shape}")
                print(f"Shape of output mask: {data_comb.shape}")

            # Ensure the shapes match before combining
            if data.shape == data_comb.shape:
                # Apply "0 wins" logic: if any mask has 0 in a pixel, the combined mask will have 0
                data_comb[(data == 0.)] = 0.
            else:
                print(f"Error: Could not combine mask {input_mask} due to mismatched shapes.")

    # Write the combined data to the output mask file
    pyfits.writeto(output_mask, data_comb, overwrite=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combine multiple masks using pixel-by-pixel 'AND' operation where 0 dominates")
    
    parser.add_argument("-o", "--output_mask", help="Name of the output mask file (default: mask.fits)", default=None)
    parser.add_argument("input_masks", nargs="+", help="Paths to the input mask files")
    parser.add_argument("-i", "--invert_masks", nargs="*", help="List of masks to invert by index (e.g., '0 2' to invert first and third mask)", default=[])

    args = parser.parse_args()

    # Convert invert_masks to boolean list
    invert_flags = [False] * len(args.input_masks)
    for idx in args.invert_masks:
        invert_flags[int(idx)] = True

    # Combine masks with the inversion option and 0 overwrites 1 behavior
    combine_masks(args.output_mask, args.input_masks, invert_flags)
