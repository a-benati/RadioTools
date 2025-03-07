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

def get_central_ra_dec(header):
    ra = header.get('CRVAL1')
    dec = header.get('CRVAL2')
    if ra is None or dec is None:
        raise ValueError('Central RA/Dec not found in the header file.')
    return ra, dec

def plot_beam_pattern(beam, freq, pol, outputplot, central_ra, central_dec):
    # Generate a grid of distances from the center
    r = np.linspace(-3, 3, 100)  # -3 to 3 degrees
    x = r * np.cos(0)
    y = r * np.sin(0)
    pattern = getattr(beam, pol)
    beam_values = pattern(x, y, freq / 1e6)
    
    # Normalize the beam pattern
    beam_values /= np.amax(beam_values)
    
    # Find the FWHM
    half_max = 0.5
    indices = np.where(np.isclose(beam_values, half_max, atol=0.1))[0]
    print("len(indices) = ", len(indices))
    if len(indices) > 0:
        fwhm_r = r[indices[0]] * 2  # FWHM is twice the first r value where beam response equals half_max
        beam_response_at_fwhm = beam_values[indices[0]]
        fwhm_position_deg = r[indices[0]]
        print(f"FWHM: {fwhm_r:.4f} degrees")
        print(f"Beam response at FWHM: {beam_response_at_fwhm:.4f}")

        # Calculate the RA and Dec at the FWHM
        fwhm_position_ra = central_ra + fwhm_position_deg / np.cos(np.radians(central_dec))
        fwhm_position_dec = central_dec + fwhm_position_deg
        print(f"FWHM position in RA: {fwhm_position_ra:.6f} degrees")
        print(f"FWHM position in Dec: {fwhm_position_dec:.6f} degrees")
    else:
        fwhm_r = None
        beam_response_at_fwhm = None
        fwhm_position_ra = None
        fwhm_position_dec = None
        print("FWHM could not be determined from the beam pattern.")
    
    plt.figure()
    plt.plot(r, beam_values, label=f'{pol}-polarization at {freq/1e6} MHz')
    plt.axhline(half_max, color='red', linestyle='--', label='Half Maximum')
    if fwhm_r:
        plt.axvline(fwhm_r / 2, color='green', linestyle='--', label=f'FWHM = {fwhm_r:.4f} deg')
        plt.axvline(-fwhm_r / 2, color='green', linestyle='--', label=f'FWHM = {-fwhm_r:.4f} deg')
    plt.xlabel('Distance from Center (degrees)')
    plt.ylabel('Beam Response')
    plt.title('MeerKAT Primary Beam Pattern')
    plt.legend(loc='lower left')
    plt.grid(True)
    plt.savefig(outputplot)
    plt.close()
    print(f"Beam pattern plot saved as {outputplot}")

def pb_correction(L, UHF, pol, savebeam, outputcorr, outputbeam, plotbeam, outputplot, *images):
    if L == True:
        beam = JimBeam('MKAT-AA-L-JIM-2020')
    elif UHF == True:
        beam = JimBeam('MKAT-AA-UHF-JIM-2020')
    else:
        raise ValueError('You must specify either L-band or UHF-band beam.')

    for image in images:
        print(f"Processing {image}...")
        with fits.open(image) as f:
            data = f[0].data
            header = f[0].header
            # Find the frequency from the header file
            freq = get_freq(header)
            print(freq)
            # Get the central RA and Dec from the header file
            central_ra, central_dec = get_central_ra_dec(header)
            print(f"Central RA: {central_ra}, Central Dec: {central_dec}")
            # Get the pixel size in degrees
            pix2deg = abs(header['CDELT1'])
            # Get the image size in degrees
            image_size = data.shape[2] * pix2deg
            # Generate a grid in degrees
            margin = np.linspace(- image_size / 2, image_size / 2, data.shape[-1])
            x, y = np.meshgrid(margin, margin)
            # Get the beam pattern
            pattern = getattr(beam, pol)
            beam_pattern = pattern(x, y, freq/1e6)
            # Normalize the beam pattern
            beam_pattern = beam_pattern / np.amax(beam_pattern)
            # Apply the beam correction
            data[:, :, :, :] = data[:, :, :, :] / beam_pattern
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
            if plotbeam:
                # Plot the beam pattern
                plot_beam_pattern(beam, freq, pol, outputplot, central_ra, central_dec)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Apply the primary beam correction to MeerKAT data either in the L or in the UHF band.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-L', dest='L', action='store_true', help='Use the L-band beam')
    group.add_argument('-UHF', dest='UHF', action='store_true', help='Use the UHF-band beam')
    parser.add_argument('--pol', dest='pol', default='I', help='Polarization that should be used. Options: "I", "H", "V" (default: "I")')
    parser.add_argument('--savebeam', dest='savebeam', action='store_true', help='Save the beam file (default: False)')
    parser.add_argument('--outputcorr', dest='outputcorr', default='pbcor', help='Name of the suffix for the corrected fits file (default: {input}_pbcor.fits)')
    parser.add_argument('--outputbeam', dest='outputbeam', default='beam', help='Name of the suffix for the output beam fits file (default: {input}_beam.fits)')
    parser.add_argument('--plotbeam', dest='plotbeam', action='store_true', help='Plot the beam pattern (default: False)')
    parser.add_argument('--outputplot', dest='outputplot', default='beam_pattern.png', help='Filename for the output beam pattern plot (default: beam_pattern.png)')
    parser.add_argument('images', nargs='+', help='List of input images.')

    args = parser.parse_args()

    pb_correction(args.L, args.UHF, args.pol, args.savebeam, args.outputcorr, args.outputbeam, args.plotbeam, args.outputplot, *args.images)
