#!/usr/bin/env python3

import argparse
import numpy as np
import os

import sys
sys.path.append('/homes/a.benati/packages/scripts/')
from lib_fits import Image

def reg2fits(regionfile, imgfile, output=None):
    im = Image(imgfile)
    im.apply_region(regionfile, blankvalue=1, invert=False)
    im.apply_region(regionfile, blankvalue=0, invert=True)

    if output is None:
        output_filename = os.path.splitext(imgfile)[0] + '_mask.fits'
    else:
        output_filename = output

    im.write(output_filename)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert a region file to a FITS mask")
    parser.add_argument("-reg", "--regionfile", required=True, help="Path to the region file")
    parser.add_argument("-img", "--imgfile", required=True, help="Path to the FITS image file")
    parser.add_argument("-o", "--output", help="Name of the output FITS mask file (default: {input}_mask.fits)", default=None)

    args = parser.parse_args()

    reg2fits(args.regionfile, args.imgfile, args.output)