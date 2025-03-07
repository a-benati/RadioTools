#!/usr/bin/env python3
import argparse
import os

def apply_to_channels(num_channels, command):
    for i in range(num_channels):
        # Format the channel index as a 4-digit number
        channel_str = f"{i:04}"
        
        # Replace '*' in the command with the current channel index
        channel_command = command.replace('*', channel_str)
        
        # Execute the command for the current channel
        print(f"Applying command to channel {channel_str}: {channel_command}")
        os.system(channel_command)

def main():
    parser = argparse.ArgumentParser(description="Apply a command to multiple channels")
    
    # Add arguments for the number of channels and the command
    parser.add_argument("-n", "--num_channels", type=int, required=True, 
                        help="Number of channels to process")
    parser.add_argument("-c", "--command", type=str, required=True, 
                        help="Command to apply to each channel, with '*' to be replaced by channel index")
    
    # Parse the arguments
    args = parser.parse_args()
    
    # Apply the command to each channel
    apply_to_channels(args.num_channels, args.command)

if __name__ == "__main__":
    main()
