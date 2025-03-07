#!/usr/bin/env python3

import argparse
import glob
from astropy.io import fits
import numpy as np

def get_beam_info(fits_file):
    with fits.open(fits_file) as hdul:
        header = hdul[0].header
        bmaj = header.get('BMAJ', np.nan) * 3600  # Convert to arcsec
        bmin = header.get('BMIN', np.nan) * 3600  # Convert to arcsec
    return bmaj, bmin

def filter_channels(beams, tolerance=3):
    bmaj_values = np.array([b[0] for b in beams])
    bmin_values = np.array([b[1] for b in beams])
    median_bmaj, median_bmin = np.median(bmaj_values), np.median(bmin_values)
    
    valid_channels = []
    invalid_channels = []
    for i, (bmaj, bmin) in enumerate(beams):
        if abs(bmaj - median_bmaj) <= tolerance and abs(bmin - median_bmin) <= tolerance:
            valid_channels.append(i)
        else:
            invalid_channels.append(i)
    
    return valid_channels, invalid_channels, max(bmaj_values[valid_channels]), max(bmin_values[valid_channels])

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('base_name', type=str, help='Base name of the FITS images')
    parser.add_argument('-n', '--num_channels', type=int, help='Number of channels')
    args = parser.parse_args()
    
    beams = []
    fits_files = []
    for i in range(args.num_channels):
        fits_file = f"{args.base_name}-{i:04d}-I-image.fits"
        fits_files.append(fits_file)
        bmaj, bmin = get_beam_info(fits_file)
        beams.append((bmaj, bmin))
    
    valid_channels, invalid_channels, max_bmaj, max_bmin = filter_channels(beams)
    
    print("All channels:")
    for i, (bmaj, bmin) in enumerate(beams):
        print(f"Channel {i:04d}: BMAJ = {bmaj:.2f} arcsec, BMIN = {bmin:.2f} arcsec")
    
    print("\nValid channels:")
    for i in valid_channels:
        print(f"Channel {i:04d}: BMAJ = {beams[i][0]:.2f} arcsec, BMIN = {beams[i][1]:.2f} arcsec")
    
    print("\nExcluded channels:")
    for i in invalid_channels:
        print(f"Channel {i:04d}: BMAJ = {beams[i][0]:.2f} arcsec, BMIN = {beams[i][1]:.2f} arcsec")

    print(f"\nMax BMAJ: {max_bmaj:.2f} arcsec")
    print(f"Max BMIN: {max_bmin:.2f} arcsec")
    
if __name__ == '__main__':
    main()
