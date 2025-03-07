#! /usr/bin/env python

import numpy as np
import casacore.tables as ct
import argparse

def calculate_uv_distances(ms_file):
    # Open the measurement set
    tb = ct.table(ms_file, ack=False)

    # Get the UVW coordinates
    uvw = tb.getcol('UVW')

    # Calculate the uv-distance in meters
    uv_distances = np.sqrt(uvw[0]**2 + uvw[1]**2)

    tb.close()
    return uv_distances

def main(ms_file_1280, ms_file_816):
    # Calculate uv-distances for both measurement sets
    uv_distances_1280 = calculate_uv_distances(ms_file_1280)
    uv_distances_816 = calculate_uv_distances(ms_file_816)

    # Convert frequencies to wavelengths
    c = 3e8  # Speed of light in m/s
    lambda_1280 = c / 1280e6
    lambda_816 = c / 816e6

    # Calculate uv-distances in wavelengths
    uv_distances_lambda_1280 = uv_distances_1280 / lambda_1280
    uv_distances_lambda_816 = uv_distances_816 / lambda_816

    # Find the common uv-range in wavelengths
    min_lambda = max(np.min(uv_distances_lambda_1280), np.min(uv_distances_lambda_816))
    max_lambda = min(np.max(uv_distances_lambda_1280), np.max(uv_distances_lambda_816))

    # Calculate the number of baselines and flagged baselines
    num_baselines = len(uv_distances_1280)
    num_flagged_baselines_1280 = len(uv_distances_1280) - np.count_nonzero(uv_distances_1280)
    num_flagged_baselines_816 = len(uv_distances_816) - np.count_nonzero(uv_distances_816)

    # Print results
    print(f"1280 MHz: Min Lambda: {np.min(uv_distances_lambda_1280):.2f}, Max Lambda: {np.max(uv_distances_lambda_1280):.2f}")
    print(f"816 MHz: Min Lambda: {np.min(uv_distances_lambda_816):.2f}, Max Lambda: {np.max(uv_distances_lambda_816):.2f}")
    print(f"Common uv-range in lambda: {min_lambda:.2f} to {max_lambda:.2f}")
    print(f"Total number of baselines: {num_baselines}")
    print(f"Number of flagged baselines at 1280 MHz: {num_flagged_baselines_1280}")
    print(f"Number of flagged baselines at 816 MHz: {num_flagged_baselines_816}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate uv-range and baseline statistics")
    parser.add_argument("ms_file_1280", help="Measurement set file for 1280 MHz")
    parser.add_argument("ms_file_816", help="Measurement set file for 816 MHz")
    args = parser.parse_args()

    main(args.ms_file_1280, args.ms_file_816)
