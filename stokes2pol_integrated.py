#!/usr/bin/env python3

import argparse
import numpy as np
from astropy.io import fits
import os
import subprocess

def compute_integrated_polarisation_fraction(I, Q, U):
    """Compute the integrated polarisation fraction as sqrt(Q_int^2 + U_int^2) / I_int."""
    Q_int = np.sum(Q)
    U_int = np.sum(U)
    I_int = np.sum(I)
    
    P_lin_int = np.sqrt(Q_int**2 + U_int**2)
    frac_pol_int = P_lin_int / I_int if I_int != 0 else 0
    return frac_pol_int

def compute_integrated_polarisation_angle(Q, U):
    """Compute the integrated polarisation angle in degrees as 0.5 * arctan(U_int / Q_int)."""
    Q_int = np.sum(Q)
    U_int = np.sum(U)
    
    pol_angle_int = 0.5 * np.arctan2(U_int, Q_int)  # Arctan2 handles the quadrant correctly
    pol_angle_deg_int = np.degrees(pol_angle_int)
    return pol_angle_deg_int

def compute_total_polarisation(Q, U):
    """Compute the total polarisation P as sqrt(Q^2 + U^2) for each pixel."""
    P = np.sqrt(Q**2 + U**2)
    return P

def compute_peak_ratio(P, I):
    """Compute the ratio between the peak P and peak I values."""
    peak_P = np.max(P)
    peak_I = np.max(I)
    ratio = peak_P / peak_I if peak_I != 0 else 0
    return ratio

def load_fits_data(file_path):
    """Load FITS file data."""
    with fits.open(file_path) as hdulist:
        data = hdulist[0].data
        header = hdulist[0].header
    return data, header

def generate_mask(region_file, img_file, output_mask):
    """Call the reg2fits.py script to convert a region file to a FITS mask."""
    reg2fits_path = "/RadioTools/reg2fits.py"
    cmd = ["python3", reg2fits_path, "-reg", region_file, "-img", img_file, "-o", output_mask]
    subprocess.run(cmd, check=True)
    print(f"Mask generated: {output_mask}")

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Compute integrated polarisation fraction, angle, P, and peak P/I ratio from Stokes I, Q, U fits files inside a given region")
    parser.add_argument("--I", required=True, help="FITS file for Stokes I")
    parser.add_argument("--Q", required=True, help="FITS file for Stokes Q")
    parser.add_argument("--U", required=True, help="FITS file for Stokes U")
    parser.add_argument("--region", required=True, help="Region file for defining the area of interest")
    parser.add_argument("--output_prefix", required=True, help="Output prefix for the results")

    args = parser.parse_args()

    # Load data from FITS files
    I_data, I_header = load_fits_data(args.I)
    Q_data, Q_header = load_fits_data(args.Q)
    U_data, U_header = load_fits_data(args.U)

    # Generate mask from the region file
    mask_file = f"{args.output_prefix}_mask.fits"
    generate_mask(args.region, args.I, mask_file)

    # Load the generated mask
    mask_data, _ = load_fits_data(mask_file)

    # Apply the mask to the Stokes parameters
    I_data_masked = I_data * mask_data
    Q_data_masked = Q_data * mask_data
    U_data_masked = U_data * mask_data

    # Compute integrated polarisation fraction, angle, total polarisation P, and peak P/I ratio
    pol_fraction_int = compute_integrated_polarisation_fraction(I_data_masked, Q_data_masked, U_data_masked)
    pol_angle_int = compute_integrated_polarisation_angle(Q_data_masked, U_data_masked)
    P_total = compute_total_polarisation(Q_data_masked, U_data_masked)
    P_total_sum = np.sum(P_total)

    # Compute the ratio between the peak P and peak I values
    peak_ratio = compute_peak_ratio(P_total, I_data_masked)

    # Output the integrated polarisation fraction, angle, P, and peak P/I ratio
    print(f"Integrated Polarisation Fraction: {pol_fraction_int}")
    print(f"Integrated Polarisation Angle: {pol_angle_int} degrees")
    print(f"Total Polarisation (P): {P_total_sum}")
    print(f"Peak P/I Ratio: {peak_ratio}")

    # Save the results to a text file
    output_file = f"{args.output_prefix}_integrated_polarisation_results.txt"
    with open(output_file, 'w') as f:
        f.write(f"Integrated Polarisation Fraction: {pol_fraction_int}\n")
        f.write(f"Integrated Polarisation Angle: {pol_angle_int} degrees\n")
        f.write(f"Total Polarisation (P): {P_total_sum}\n")
        f.write(f"Peak P/I Ratio: {peak_ratio}\n")
    
    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    main()
