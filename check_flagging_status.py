#!/usr/bin/env python

import argparse
from casacore.tables import table

def check_flagging_status(ms_path):
    # Open the MS file
    ms = table(ms_path, readonly=True)

    # Get the flagging columns
    flag_cols = ms.colnames()

    # Iterate over flagging columns and print flagging statistics
    for col in flag_cols:
        flags = ms.getcol(col)
        num_flags = flags.sum()
        total_data = flags.size
        percent_flagged = (num_flags / total_data) * 100
        print(f"Column: {col}")
        print(f"Number of flagged data points: {num_flags}")
        print(f"Total data points: {total_data}")
        print(f"Percentage flagged: {percent_flagged:.2f}%")
        print()

    # Close the MS file
    ms.close()

def main():
    # Create an argument parser
    parser = argparse.ArgumentParser(description="Check flagging status of a Measurement Set (MS) file")

    # Add the argument for the path to the MS file
    parser.add_argument("ms_file", help="Path to the MS file to check flagging status")

    # Parse the command-line arguments
    args = parser.parse_args()

    # Call the function to check flagging status
    check_flagging_status(args.ms_file)

if __name__ == "__main__":
    main()
