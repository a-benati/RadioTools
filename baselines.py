#! /usr/bin/env python

import numpy as np
from casacore.tables import table
import argparse

def calculate_baselines(ms_file, uvcut_lambda):
    # Constants
    frequency = 816e6  # Frequency in Hz
    c = 3e8  # Speed of light in m/s

    # Calculate the wavelength
    wavelength = c / frequency

    # Open the ANTENNA table
    antenna_table = table(ms_file + '/ANTENNA')
    positions = antenna_table.getcol('POSITION')
    antenna_names = antenna_table.getcol('NAME')
    antenna_table.close()

    # Calculate the baseline lengths
    def baseline_length(pos1, pos2):
        return np.sqrt(np.sum((pos1 - pos2) ** 2))

    baselines = []
    num_antennas = positions.shape[0]
    for i in range(num_antennas):
        for j in range(i + 1, num_antennas):
            length = baseline_length(positions[i], positions[j])
            baselines.append((length, antenna_names[i], antenna_names[j]))

    baselines = np.array(baselines, dtype=[('length', 'f8'), ('antenna1', 'U50'), ('antenna2', 'U50')])

    # Convert baseline lengths to wavelengths
    baselines_lambda = baselines['length'] / wavelength

    # Filter baselines based on uvcut_lambda
    cut_baselines = baselines_lambda[baselines_lambda < uvcut_lambda]
    kept_baselines = baselines_lambda[baselines_lambda >= uvcut_lambda]

    # Calculate statistics
    min_baseline = np.min(baselines['length'])
    max_baseline = np.max(baselines['length'])

    min_baseline_info = baselines[np.argmin(baselines['length'])]
    max_baseline_info = baselines[np.argmax(baselines['length'])]

    # Open the MS table to check flagged data
    ms_table = table(ms_file, ack=False)
    flagged_data = ms_table.getcol('FLAG')
    antenna1 = ms_table.getcol('ANTENNA1')
    antenna2 = ms_table.getcol('ANTENNA2')
    ms_table.close()

    # Calculate the number of baselines
    total_baselines = len(baselines)

    # Check for flagged baselines
    num_time_samples = flagged_data.shape[0]
    flagged_baselines = 0

    unique_baselines = set(zip(antenna1, antenna2))
    for ant1, ant2 in unique_baselines:
        baseline_indices = np.where((antenna1 == ant1) & (antenna2 == ant2))
        if np.all(flagged_data[baseline_indices]):
            flagged_baselines += 1

    unflagged_baselines = total_baselines - flagged_baselines

    print(f"Minimum baseline length: {min_baseline:.2f} meters (between {min_baseline_info['antenna1']} and {min_baseline_info['antenna2']})")
    print(f"Maximum baseline length: {max_baseline:.2f} meters (between {max_baseline_info['antenna1']} and {max_baseline_info['antenna2']})")
    print(f"Total number of unique baselines: {total_baselines}")
    print(f"Number of fully flagged baselines: {flagged_baselines}")
    print(f"Number of unflagged baselines: {unflagged_baselines}")
    print(f"Number of baselines shorter than {uvcut_lambda} lambda: {len(cut_baselines)}")
    print(f"Number of baselines longer than or equal to {uvcut_lambda} lambda: {len(kept_baselines)}")

    # Print cut baselines details
    for baseline in baselines[baselines_lambda < uvcut_lambda]:
        print(f"Cut baseline: {baseline['length']:.2f} meters (between {baseline['antenna1']} and {baseline['antenna2']})")

    return min_baseline, max_baseline, total_baselines, flagged_baselines, unflagged_baselines

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate baseline statistics from a Measurement Set (MS) file.")
    parser.add_argument("ms_file", help="Path to the Measurement Set (MS) file.")
    parser.add_argument("uvcut_lambda", type=float, help="UV cut-off in wavelengths.")
    args = parser.parse_args()

    calculate_baselines(args.ms_file, args.uvcut_lambda)
