#!/usr/bin/env python3

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
import pyregion

def apply_region_mask(data, region_file, header):
    """
    Apply a DS9 region mask to the data using pyregion.
    """
    regions = pyregion.open(region_file)
    
    # Create the mask from the region
    mask = regions.get_mask(header=header, shape=data.shape[-2:])
    
    return mask

def extract_region_data(data, mask):
    """
    Extract data using a region mask.
    """
    # Apply the mask to the data and return the region data
    masked_data = np.ma.masked_array(data, ~mask)
    return masked_data

def derive_polarisation_values(pol_fraction_data, pol_angle_data, mask):
    """
    Derive the polarisation fraction and angle for the source using the region mask.
    """
    # Apply mask to the data
    pol_fraction_region = extract_region_data(pol_fraction_data[0, 0, :, :], mask)
    pol_angle_region = extract_region_data(pol_angle_data[0, 0, :, :], mask)

    # Debug: Check if the region is empty
    print(f"Masked data shape (Pol. Fraction): {pol_fraction_region.shape}")
    print(f"Masked data shape (Pol. Angle): {pol_angle_region.shape}")

    if np.ma.is_masked(pol_fraction_region) and pol_fraction_region.count() == 0:
        print("Warning: Extracted region has no valid data for polarisation fraction.")
        return np.nan, np.nan

    if np.ma.is_masked(pol_angle_region) and pol_angle_region.count() == 0:
        print("Warning: Extracted region has no valid data for polarisation angle.")
        return np.nan, np.nan

    # Compute the median (or mean) values within the masked region
    pol_fraction = np.ma.median(pol_fraction_region)
    std_pol_fraction = np.std(pol_fraction_data, ddof=1)
    pol_angle = np.ma.median(pol_angle_region)
    std_pol_angle = np.std(pol_angle_data, ddof=1)
    
    return pol_fraction, pol_angle, std_pol_fraction, std_pol_angle

def main(polarisation_fraction_fits, polarisation_angle_fits, region_file):
    # Load the polarisation fraction and angle data
    pol_fraction_hdulist = fits.open(polarisation_fraction_fits)
    pol_angle_hdulist = fits.open(polarisation_angle_fits)
    
    pol_fraction_data = pol_fraction_hdulist[0].data
    pol_angle_data = pol_angle_hdulist[0].data
    
    # Load the WCS information from the FITS header
    wcs = WCS(pol_fraction_hdulist[0].header)
    
    # Apply region mask
    mask = apply_region_mask(pol_fraction_data, region_file, pol_fraction_hdulist[0].header)
    
    # Derive the polarisation fraction and angle at the source position
    pol_fraction, pol_angle, std_pol_fraction, std_pol_angle = derive_polarisation_values(pol_fraction_data, pol_angle_data, mask)
    
    if np.isnan(pol_fraction) or np.isnan(pol_angle):
        print("Polarisation Fraction or Angle could not be determined due to lack of valid data.")
    else:
        print(f"Polarisation Fraction: {pol_fraction * 100:.2f} ± {std_pol_fraction:.2f}%")
        print(f"Polarisation Angle: {pol_angle:.2f} ± {std_pol_angle:.2f} degrees")
    
    pol_fraction_hdulist.close()
    pol_angle_hdulist.close()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Polarisation analysis tool using DS9 regions")
    parser.add_argument("pol_angle_fits", help="FITS file with polarisation angle image")
    parser.add_argument("pol_frac_fits", help="FITS file with polarisation fraction image")
    parser.add_argument("--region", type=str, required=True, help="DS9 region file for the source region")
    
    args = parser.parse_args()
    main(args.pol_frac_fits, args.pol_angle_fits, args.region)
