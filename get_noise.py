#!/usr/bin/env python3

import argparse
import numpy as np
from astropy.io import fits
from photutils.background import MADStdBackgroundRMS

def estimate_noise(image_file):
    """Stima il rumore di un'immagine FITS usando il metodo MAD."""
    # Apri il file FITS
    hdu = fits.open(image_file)
    data = hdu[0].data
    hdu.close()

    # Converti NaN e infiniti in zero per evitare errori
    data = np.nan_to_num(data, nan=0.0, posinf=0.0, neginf=0.0)
    print(data.shape)

    # Calcola il noise con MAD (Median Absolute Deviation)
    bkg_estimator = MADStdBackgroundRMS()

    # Check if the FITS file is a cube (3D) or a single image (2D)
    if data.ndim == 2:
        noise = bkg_estimator(data)
        print(f"Noise of {image_file}: {noise:.6f} Jy/beam")
    elif data.ndim == 3:
        noise_values = np.array([bkg_estimator(slice_) for slice_ in data])
        overall_noise = np.median(noise_values)  # Median of MADs
        print(f"Noise values per slice: {noise_values}")
        print(f"Overall noise (median of MADs) of {image_file}: {overall_noise:.6f} Jy/beam")
    else:
        print("Unsupported FITS format: expected 2D or 3D data.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Estimate noise in a FITS image using MAD.")
    parser.add_argument("image", type=str, help="Path to the FITS file")
    args = parser.parse_args()

    estimate_noise(args.image)
