#!/usr/bin/env python3

import argparse
from astropy.io import fits

# Configura argparse per accettare input dalla riga di comando
parser = argparse.ArgumentParser(description="Print the header file of a FITS image")
parser.add_argument("imagename", type=str, help="Path of the FITS image")
args = parser.parse_args()

# Leggi il nome dell'immagine dalla riga di comando
imagename = args.imagename

# Apri e leggi l'immagine FITS
image = fits.open(imagename)

# Stampa l'header dell'immagine FITS
print(f"Header of the image {imagename}:")
print(image[0].header)
