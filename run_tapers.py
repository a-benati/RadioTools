#!/usr/bin/env python3

import argparse
import os
import subprocess

def run_wsclean(command, tapers):
    # Split the command into a list to work with it
    command_list = command.split()

    # Identify the index of the name parameter and remove it from the list temporarily
    if '-name' in command_list:
        name_index = command_list.index('-name') + 1
        base_name = command_list[name_index]
    else:
        print("Error: '-name' parameter is missing from the command.")
        return

    # Get the base directory and filename from the name
    base_dir = os.path.dirname(base_name)
    base_filename = os.path.basename(base_name)

    # Loop over each taper value and modify the command accordingly
    for taper in tapers:
        print(f"Running WSClean with taper: {taper}...")

        # Create the new directory name based on the taper value
        taper_dir = os.path.join(base_dir, f"taper_{taper}")
        if not os.path.exists(taper_dir):
            os.makedirs(taper_dir)

        # Create the new name by appending "_t{taper}" to the filename and putting it in the taper directory
        new_name = os.path.join(taper_dir, f"{base_filename}_t{taper}")
        command_list[name_index] = new_name

        # Add the taper option "-taper-gaussian X" to the command
        if '-taper-gaussian' in command_list:
            taper_index = command_list.index('-taper-gaussian') + 1
            command_list[taper_index] = str(taper)  # Replace existing taper value
        else:
            command_list.insert(1, '-taper-gaussian')  # Insert the taper option
            command_list.insert(2, str(taper))

        # Run the modified command and wait for it to finish
        try:
            result = subprocess.run(command_list, check=True)
            if result.returncode != 0:
                print(f"Command failed with return code: {result.returncode}")
        except subprocess.CalledProcessError as e:
            print(f"Command failed with error: {e}")
            continue  # Continue with the next taper value

def main():
    parser = argparse.ArgumentParser(description="Run WSClean with different tapers.")
    parser.add_argument('-c', '--command', type=str, required=True, help='The base WSClean command')
    parser.add_argument('-t', '--tapers', type=int, nargs='+', required=True, help='List of taper values (e.g., 25 50 75)')

    args = parser.parse_args()

    # Run the WSClean command with the specified tapers
    run_wsclean(args.command, args.tapers)

if __name__ == '__main__':
    main()
