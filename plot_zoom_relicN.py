#!/usr/bin/env python3

import sys, os 
from astropy.wcs import WCS
from astropy.visualization import MinMaxInterval, SqrtStretch, ImageNormalize, LogStretch
from astropy.io import fits
import numpy as np
from astropy import units as u
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt

sys.path.append('/packages/scripts/')
from lib_fits import flatten
import lib_plot

# Common variables
center = [-57.389218710,-56.652536022]
size = [0.385657*2,0.217885]
z = 0.0553  # redshift
kpc = 1000  # scalebar length in kpc

# L band image
def plot_image(image, vmax, ticks, output_file):
    header, radio = flatten(image)
    wcs = WCS(header)
    radio *= 1e3
    
    fig = plt.figure(figsize=(20, 21))
    ax = fig.add_subplot(1, 1, 1, projection=wcs, slices=('x', 'y'))
    lon = ax.coords['ra']
    lat = ax.coords['dec']
    lon.set_axislabel('Right Ascension (J2000)', size=20)
    lat.set_axislabel('Declination (J2000)', size=20)
    lon.set_major_formatter('hh:mm:ss')
    lat.set_major_formatter('dd:mm')
    lat.set_ticklabel(rotation=90)  # to turn dec vertical

    # Set limits
    xlims = [center[0] - size[0] / 2., center[0] + size[0] / 2.]
    ylims = [center[1] - size[1] / 2., center[1] + size[1] / 2.]
    x, y = wcs.wcs_world2pix(xlims*u.deg, ylims*u.deg, 1, ra_dec_order=True)
    ax.set_xlim(x[1], x[0])
    ax.set_ylim(y[0], y[1])
    lon.set_ticklabel(size=20)
    lat.set_ticklabel(size=20)

    # Add scalebar
    lib_plot.addScalebar(ax, wcs, z, kpc, fontsize=20, color='white')

    # Normalization and image plotting
    norm = ImageNormalize(radio, vmin=0.003, vmax=vmax, stretch=LogStretch())
    im = ax.imshow(radio, cmap='inferno', origin='lower', norm=norm)

    # Colorbar
    cbaxes = fig.add_axes([0.127, 0.72, 0.772, 0.02])  # Adjusted for proper colorbar positioning
    cbar = fig.colorbar(im, cax=cbaxes, ticks=ticks, format='$%.2f$', orientation='horizontal')
    cbar.set_label(r'Surface brightness (mJy beam$^{-1}$)', labelpad=20)
    cbaxes.xaxis.tick_top()
    cbaxes.tick_params(labelsize=20)
    cbaxes.xaxis.label.set_size(20)
    cbaxes.xaxis.set_label_position('top')

    # Save figure
    fig.savefig(output_file, bbox_inches='tight')
    plt.close(fig)

# Plot L band
# plot_image(
#     image='/data/abell_3667/l_band/node34-MFS-image-PBcorr.fits',
#     vmax=1.0,
#     ticks=[0.01, 0.1, 0.3, 0.8],
#     output_file='/data/abell_3667/img_pdf/abell3667_L_zoom_relicN.pdf'
# )

# Plot UHF band
plot_image(
    image='/data/abell_3667/img_uhf/ddcal_005-MFS-image_pbcor.fits',
    vmax=3.0,
    ticks=[0.03, 0.3, 1.0, 2.5],
    output_file='/data/abell_3667/img_pdf/abell3667_UHF_zoom_relicN.pdf'
)
