#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS

# Load the radio image data (corrected and uncorrected)
hdu_corrected = fits.open('/data/abell_3667/ddcal_pipeline/extract/facetselfcal_ddcal_extract_002-MFS-image_pbcor.fits')
#hdu_corrected = fits.open('/data/abell_3667/ddcal_pipeline/extract/facetselfcal_ddcal_extract_002-MFS-image-pb.fits')
data_corrected = hdu_corrected[0].data
wcs = WCS(hdu_corrected[0].header)

hdu_uncorrected = fits.open('/data/abell_3667/ddcal_pipeline/extract/facetselfcal_ddcal_extract_002-MFS-image.fits')
data_uncorrected = hdu_uncorrected[0].data

# Check the shape of the data
print(f"Shape of corrected image data: {data_corrected.shape}")
print(f"Shape of uncorrected image data: {data_uncorrected.shape}")

# Extract a 2D slice if necessary
if len(data_corrected.shape) > 2:
    # Assuming we want the first slice if there are more than 2 dimensions
    data_corrected = data_corrected[0, 0]
    data_uncorrected = data_uncorrected[0, 0]
    print(f"Using first slice of the data. New shape: {data_corrected.shape}")

# Define the central frequency and calculate the FWHM of the primary beam
central_frequency_mhz = 816
dish_diameter_m = 13.5

def calculate_fwhm(frequency_mhz, dish_diameter_m):
    c = 3e8
    frequency_hz = frequency_mhz * 1e6
    wavelength_m = c / frequency_hz
    k = 1.2
    fwhm_radians = k * wavelength_m / dish_diameter_m
    fwhm_degrees = np.degrees(fwhm_radians)
    return fwhm_degrees

fwhm_degrees = -1.7477#calculate_fwhm(central_frequency_mhz, dish_diameter_m)
half_fwhm_degrees = fwhm_degrees / 2

print(f"FWHM: {fwhm_degrees:.2f} degrees")
print(f"Half FWHM: {half_fwhm_degrees:.2f} degrees")

# Calculate the pixel distance corresponding to half FWHM in degrees
# Assuming the image is centered and square
pixel_scale = np.abs(wcs.wcs.cdelt[0])  # degrees per pixel
half_fwhm_pixels = half_fwhm_degrees / pixel_scale

# Center of the image
center_x, center_y = np.array(data_corrected.shape) // 2

# Create a mask for the region at half FWHM
y, x = np.ogrid[:data_corrected.shape[0], :data_corrected.shape[1]]
mask_at_half_fwhm = np.isclose((x - center_x)**2 + (y - center_y)**2, half_fwhm_pixels**2, atol=1)

# Calculate the average flux density at half FWHM for both images
avg_flux_corrected_at_half = np.mean(data_corrected[mask_at_half_fwhm])
avg_flux_uncorrected_at_half = np.mean(data_uncorrected[mask_at_half_fwhm])

print(f"Average flux density at half FWHM in corrected image: {avg_flux_corrected_at_half}")
print(f"Average flux density at half FWHM in uncorrected image: {avg_flux_uncorrected_at_half}")

# Check if the corrected flux density is approximately double at half FWHM
expected_flux_increase_at_half = 2 * avg_flux_uncorrected_at_half
is_corrected_flux_appropriate_at_half = np.isclose(avg_flux_corrected_at_half, expected_flux_increase_at_half, rtol=0.1)

print(f"Is the corrected flux density at half FWHM approximately double? {is_corrected_flux_appropriate_at_half}")

if is_corrected_flux_appropriate_at_half:
    print("Primary beam correction seems correct at half FWHM.")
else:
    print("Primary beam correction might be incorrect at half FWHM.")