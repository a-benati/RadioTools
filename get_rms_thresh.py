#!/usr/bin/env python3

import argparse
import numpy as np
import os
from astropy.io import fits
from photutils.background import MADStdBackgroundRMS

def estimate_noise_from_residuals(base_path, num_channels):
    """Estimate the RMS noise from the residual FITS and determine the RMS threshold."""
    
    noise_per_channel = []
    channel_ids = []
    filenames = []
    
    # Estimatore del background noise
    bkg_estimator = MADStdBackgroundRMS()
    
    for i in range(num_channels):
        filename = f"{base_path}-{i:04d}-I-residual.fits"
        filenames.append(filename)
        
        if not os.path.exists(filename):
            print(f"Warning: File {filename} not found, skipping.")
            continue
        
        with fits.open(filename) as hdu:
            data = hdu[0].data
        
        # Assumiamo che i dati siano 2D (Y, X)
        channel_data = np.nan_to_num(data, nan=0.0, posinf=0.0, neginf=0.0)
        noise = bkg_estimator(channel_data)
        
        noise_per_channel.append(noise)
        channel_ids.append(i)
    
    noise_per_channel = np.array(noise_per_channel)

    if len(noise_per_channel) == 0:
        print("Error: No valid file found.")
        return
    
    # Stima la threshold ottimale
    median_noise = np.median(noise_per_channel)
    mean_noise = np.mean(noise_per_channel)
    std_noise = np.std(noise_per_channel)

    threshold_2sigma = 2 * median_noise
    threshold_90percentile = np.percentile(noise_per_channel, 90)
    threshold_3sigma = mean_noise + 3 * std_noise

    # Scegli il metodo migliore
    if std_noise / mean_noise < 0.2:
        rms_threshold = threshold_2sigma
        method_used = "2σ over the median"
    elif std_noise / mean_noise > 0.5:
        rms_threshold = threshold_3sigma
        method_used = "3σ clipping"
    else:
        rms_threshold = threshold_90percentile
        method_used = "90° percentile"

    # Identifica i canali da escludere
    excluded_channels = [channel_ids[i] for i in range(len(noise_per_channel)) if noise_per_channel[i] > rms_threshold]

    # Stampa i risultati
    print("\Channel | RMS Noise (Jy/beam)")
    print("-" * 28)
    for i, noise in zip(channel_ids, noise_per_channel):
        mark = "  <-- niose > threshold" if noise > rms_threshold else ""
        print(f"{i:04d}   | {noise:.3e}{mark}")

    print("-" * 28)
    print(f"Method: {method_used}")
    print(f"Threshold RMS: {rms_threshold:.3e} Jy/beam")
    print(f"Excluded channels: {excluded_channels}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Estimate the RMS noise from the residual FITS and determine the RMS threshold.")
    parser.add_argument("base_path", type=str, help="Path of the residual files (without the number of channel).")
    parser.add_argument("-n", "--num_channels", type=int, required=True, help="Total number of channels to be analyzed.")
    args = parser.parse_args()

    estimate_noise_from_residuals(args.base_path, args.num_channels)
