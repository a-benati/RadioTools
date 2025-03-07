#!/usr/bin/env python3

import sys
import os
from astropy.wcs import WCS
from astropy.visualization import MinMaxInterval, SqrtStretch, ImageNormalize, LinearStretch
from astropy.io import fits
import numpy as np
from astropy import units as u
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

sys.path.append('/packages/scripts/')
from lib_fits import flatten
import lib_plot

# File paths for the two relics
image_north = '/data/abell_3667/spectral_index_map/spidxmap_relicN.fits'
image_south = '/data/abell_3667/spectral_index_map/spidxmap_relicS.fits'

# Load and process the north image
header_north, radio_north = flatten(image_north)
wcs_north = WCS(header_north)

# Load and process the south image
header_south, radio_south = flatten(image_south)
wcs_south = WCS(header_south)

# Create a figure with two subplots side by side
fig = plt.figure(figsize=(30, 22))
gs = GridSpec(1, 2, width_ratios=[10, 20], wspace=0.1)

# Function to plot the images
def plot_image(ax, wcs, radio, center, size, z, kpc, norm, index):
    #ax = fig.add_subplot(gs[index], projection=wcs, slices=('x', 'y'))
    lon = ax.coords['ra']
    lat = ax.coords['dec']
    lon.set_axislabel('Right Ascension (J2000)', size=20)
    lat.set_axislabel('Declination (J2000)', size=20)
    lon.set_major_formatter('hh:mm:ss')
    lat.set_major_formatter('dd:mm')
    lat.set_ticklabel(rotation=90)

    xlims = [center[0] - size[0] / 2., center[0] + size[0] / 2.]
    ylims = [center[1] - size[1] / 2., center[1] + size[1] / 2.]
    x, y = wcs.wcs_world2pix(xlims * u.deg, ylims * u.deg, 1, ra_dec_order=True)
    ax.set_xlim(x[1], x[0])
    ax.set_ylim(y[0], y[1])
    lon.set_ticklabel(size=20)
    lat.set_ticklabel(size=20)

    lib_plot.addScalebar(ax, wcs, z, kpc, fontsize=20, color='black')

    im = ax.imshow(radio, cmap='Spectral_r', origin='lower', norm=norm)

    return im

# Plot the north image
center_north = [302.82186, -56.4840023]
size_north = [0.57 / np.cos(center_north[1] * np.pi / 180), 0.53]
norm_north = ImageNormalize(radio_north, vmin=-3., vmax=-0.8, stretch=LinearStretch())
ax_north = fig.add_subplot(gs[1], projection=wcs_north, slices=('x', 'y'))
#ax_north = fig.add_axes([0.05, 0.1, 10/30, 0.8], projection=wcs_north, slices=('x', 'y'))
im_north = plot_image(ax_north, wcs_north, radio_north, center_north, size_north, 0.0553, 1000, norm_north, 1)

# Plot the south image
center_south = [303.5453, -57.05]
size_south = [0.5, 0.53]
norm_south = ImageNormalize(radio_south, vmin=-3., vmax=-0.8, stretch=LinearStretch())
ax_south = fig.add_subplot(gs[0], projection=wcs_south, slices=('x', 'y'))
im_south = plot_image(ax_south, wcs_south, radio_south, center_south, size_south, 0.0553, 500, norm_south, 0)

# Add a combined colorbar
cbaxes = fig.add_axes([0.127, 0.83, 0.772, 0.02])
cbar = fig.colorbar(im_north, cax=cbaxes, format='$%.1f$', orientation='horizontal')
cbar.set_label(r'Spectral index', labelpad=20)
cbaxes.xaxis.tick_top()
cbaxes.tick_params(labelsize=20)
cbaxes.xaxis.label.set_size(25)
cbaxes.xaxis.set_label_position('top')

# Save the figure to a PDF file
fig.savefig('spectral_index_map.pdf', bbox_inches='tight')
