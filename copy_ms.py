#!/usr/bin/env python

import casacore.tables as pt
import os
import argparse

# Function to create a copy of the MS
def copy_ms_file(original_ms, copy_ms):
    # Check if the original MS exists
    if not os.path.exists(original_ms):
        print(f"Error: The original MS file {original_ms} does not exist.")
        return
    
    # Check if the copy already exists
    if os.path.exists(copy_ms):
        print(f"Error: The copy MS file {copy_ms} already exists. Choose a different name or delete the existing file.")
        return

    print(f"Copying {original_ms} to {copy_ms}...")
    
    # Open the original MS in read-only mode and copy it
    with pt.table(original_ms, readonly=True) as tb:
        tb.copy(copy_ms, deep=True)
    
    print(f"Copy successful! New MS file is saved as {copy_ms}")

if __name__ == "__main__":
    # Set up argparse
    parser = argparse.ArgumentParser(description="Create a deep copy of a Measurement Set (MS) file.")
    parser.add_argument("original_ms", help="Path to the original MS file")
    parser.add_argument("copy_ms", nargs='?', help="Path to save the copied MS file (default: original filename with '_copy.ms')")

    # Parse the arguments
    args = parser.parse_args()

    # If no copy name is provided, create one based on the original
    if args.copy_ms is None:
        base_name, ext = os.path.splitext(args.original_ms)
        args.copy_ms = f"{base_name}_copy{ext}"

    # Call the function to copy the MS
    copy_ms_file(args.original_ms, args.copy_ms)
