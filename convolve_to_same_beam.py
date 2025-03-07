#!/usr/bin/env python3

#!/usr/bin/env python3

import argparse
import logging
from pathlib import Path
import sys

# Add the directory containing lib_fits.py to the system path
sys.path.append('/packages/scripts')
from lib_fits import AllImages

def main():
    # Set up logging
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

    # Argument parser
    parser = argparse.ArgumentParser(description="Convolve an input FITS image to a specified beam.")
    parser.add_argument("input_image", type=str, help="Path to the input FITS image")
    
    # Beam parameters
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--beam", type=float, help="Target beam size in arcseconds (circular beam)")
    group.add_argument("--bmaj_bmin", nargs=2, type=float, metavar=('BMAJ', 'BMIN'),
                       help="Target major and minor axis beam size in arcseconds (non-circular beam)")
    
    parser.add_argument("--BPA", type=float, help="Beam position angle (degrees)", default=0.0)

    # Parse arguments
    args = parser.parse_args()

    # Verify the input image exists
    input_path = Path(args.input_image)
    if not input_path.is_file():
        logging.error(f"Input image {input_path} does not exist.")
        sys.exit(1)

    # Determine the target beam
    if args.beam:
        # Circular beam
        target_beam = [args.beam, args.beam, args.BPA]
        logging.info(f"Using a circular beam with size {args.beam} arcsec and BPA {args.BPA} degrees.")
    else:
        # Non-circular beam
        bmaj, bmin = args.bmaj_bmin
        target_beam = [bmaj, bmin, args.BPA]
        logging.info(f"Using a non-circular beam with BMAJ {bmaj} arcsec, BMIN {bmin} arcsec, and BPA {args.BPA} degrees.")

    # Convolution process
    try:
        # Initialize AllImages with the input image
        all_images = AllImages([str(input_path)])

        # Convolve to the specified beam
        all_images.convolve_to(beam=target_beam)
        all_images.write('conv', inflate=True)
        logging.info("Convolution completed successfully.")
    except Exception as e:
        logging.error(f"Error during convolution: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
