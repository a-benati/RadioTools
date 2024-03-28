#!/usr/bin/env python

from astropy.io import fits
from astropy.wcs import WCS
import argparse
import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append('/homes/a.benati/packages/katbeam/')
from katbeam import JimBeam

def get_freq(header):
    if 'FREQ' in header['CTYPE3']:
        freq = header['CRVAL3']
    elif 'FREQ' in header['CTYPE4']:
        freq = header['CRVAL4']
    else:
        raise ValueError('Frequency not found in the header file.')
    return freq

def pb_correction(L, UHF, pol, savebeam, outputcorr, outputbeam, *images):
    if L == True:
        beam = JimBeam('MKAT-AA-L-JIM-2020')

    if UHF == True:
        beam = JimBeam('MKAT-AA-UHF-JIM-2020')

    for image in images:
        print(f"Processing {image}...")
        with fits.open(image) as f:
            data = f[0].data
            header = f[0].header
            # Find the frequency from the header file
            freq = get_freq(header)
            print(freq)
            # Get the pixel size in degrees
            pix2deg = abs(header['CDELT1'])
            # Get the image size in degrees
            image_size = data.shape[2] * pix2deg
            # Generate a grid in degrees
            margin = np.linspace(- image_size / 2, image_size / 2, 10000)
            x, y = np.meshgrid(margin, margin)
            # Get the beam pattern
            pattern = getattr(beam, pol)
            beam_pattern = pattern(x, y, freq/1e6)
            # Normalize the beam pattern
            beam_pattern = beam_pattern / np.amax(beam_pattern)
            # Apply the beam correction
            data = data / beam_pattern
            f[0].data = data
            # Save the corrected file
            corr_output_file = image.replace('.fits', f'_{outputcorr}.fits')
            print(f"Saved file: {corr_output_file}")
            fits.writeto(corr_output_file, data.astype(np.float32), header, overwrite=True)
            if savebeam:
                # Save the beam file
                fits_output = image.replace('.fits',f'_{outputbeam}.fits')
                print('Saved file: %s' % fits_output)
                fits.writeto(fits_output, beam_pattern.astype(np.float32), header=header, overwrite=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Apply the primary beam correction to MeerKAT data either in the L or in the UHF band.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-L', dest='L', action='store_true', help='Use the L-band beam')
    group.add_argument('-UHF', dest='UHF', action='store_true', help='Use the UHF-band beam')
    parser.add_argument('--pol', dest='pol', default='I', help='Polarization that should be used. Options: "I", "H", "V" (default: "I")')
    parser.add_argument('--savebeam', dest='savebeam', action='store_true', help='Save the beam file (default: False)')
    parser.add_argument('--outputcorr', dest='outputcorr', default='pbcor', help='Name of the suffix for the corrected fits file (default: {input}_pbcor.fits)')
    parser.add_argument('--outputbeam', dest='outputbeam', default='beam', help='Name of ths suffix for the output beam fits file (default: {input}_beam.fits)')
    parser.add_argument('images', nargs='+', help='List of input images.')

    args = parser.parse_args()

    pb_correction(args.L, args.UHF, args.pol, args.savebeam, args.outputcorr, args.outputbeam, *args.images)