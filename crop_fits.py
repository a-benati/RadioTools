#!/usr/bin/env python3

import argparse
from astropy.io import fits
from astropy.wcs import WCS
from regions import Regions
import numpy as np

def crop_fits(fits_file, reg_file):
    # Read the FITS file
    with fits.open(fits_file) as hdul:
        hdu = hdul[0]
        wcs = WCS(hdu.header).celestial  # Get the WCS
        data = hdu.data.squeeze()  # Remove extra dimensions

    # Read the region file
    region = Regions.read(reg_file, format='ds9')[0]  # Assume only one region is present
    
    # Convert region from celestial coordinates to pixel
    region_pixel = region.to_pixel(wcs)
    
    # Obtain the bounding box of the region in pixel
    bbox = region_pixel.bounding_box
    xmin, xmax = int(bbox.ixmin), int(bbox.ixmax)
    ymin, ymax = int(bbox.iymin), int(bbox.iymax)

    # Perform the crop
    cropped_data = data[ymin:ymax, xmin:xmax]

    # Create the new header with updated WCS
    new_header = wcs.to_header()
    
    # Write the new FTIS file
    output_file = fits_file.replace(".fits", "-crop.fits")
    fits.writeto(output_file, cropped_data, new_header, overwrite=True)
    print(f"Cropped image saved as: {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Crop a FITS image using a DS9 region file.")
    parser.add_argument("fits_file", help="Input FITS file")
    parser.add_argument("reg_file", help="DS9 region file")
    args = parser.parse_args()
    
    crop_fits(args.fits_file, args.reg_file)
