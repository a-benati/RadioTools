#!/usr/bin/env python3

import numpy as np
from astropy.io import fits

# Load image
#image_file = "/data/abell_3667/imaging/sub_sources/uvcut/uvcut12000.0/diffuse_sub-MFS-residual_pbcor.fits"
image_file = "/data/abell_3667/imaging/sub_sources/uvcut/uvcut12000.0/taper_50asec/diffuse_sub_t50asec-MFS-residual_pbcor.fits"
#image_file = "/data/abell_3667/calibration/polcal/A3667/RM_SYNTH_A3667/RM_SYNTH_A3667-MFS-Q-image_pbcor.fits"
data = fits.getdata(image_file)
data = data[0, 0, :, :]

# Define region size (e.g., 200x200 pixels)
region_size = 200
half_size = region_size // 2

# Function to measure noise in a given region
def measure_noise(region):
    mean = np.mean(region)
    std_dev = np.std(region)
    return mean, std_dev

# Identify blank regions manually or define positions
positions = [(6295, 5243), (5769, 6094), (4581, 5938), (4646,5337)]#, (4854, 5369)]  # Example positions (x, y)
noise_values = []

for x, y in positions:
    region = data[y-half_size:y+half_size, x-half_size:x+half_size]
    mean, std_dev = measure_noise(region)
    noise_values.append(std_dev)

# Average noise
average_noise = np.mean(noise_values)
print(f"Estimated Noise: {average_noise} Jy/beam")
