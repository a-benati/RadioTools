#!/usr/bin/env python3

import sys
import os
from astropy.wcs import WCS
from astropy.visualization import MinMaxInterval, SqrtStretch, ImageNormalize, LogStretch
from astropy.io import fits
import numpy as np
from astropy import units as u
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt

sys.path.append('/homes/a.benati/packages/scripts/')
from lib_fits import flatten
import lib_plot

def plot_relic(image_path, center, size, scalebar_kpc, norm_params, ticks, figsize, output_file):
    header, radio = flatten(image_path)
    wcs = WCS(header)
    radio *= 1e3

    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(1, 1, 1, projection=wcs, slices=('x', 'y'))
    lon = ax.coords['ra']
    lat = ax.coords['dec']
    lon.set_axislabel('Right Ascension (J2000)', size=20)
    lat.set_axislabel('Declination (J2000)', size=20)
    lon.set_major_formatter('hh:mm:ss')
    lat.set_major_formatter('dd:mm')
    lat.set_ticklabel(rotation=90)  # to turn dec vertical

    size[0] /= np.cos(center[1] * np.pi / 180)
    xlims = [center[0] - size[0] / 2., center[0] + size[0] / 2.]
    ylims = [center[1] - size[1] / 2., center[1] + size[1] / 2.]
    x, y = wcs.wcs_world2pix(xlims * u.deg, ylims * u.deg, 1, ra_dec_order=True)
    ax.set_xlim(x[1], x[0])
    ax.set_ylim(y[0], y[1])
    lon.set_ticklabel(size=20)
    lat.set_ticklabel(size=20)

    # scalebar
    z = 0.0553  # redshift
    lib_plot.addScalebar(ax, wcs, z, scalebar_kpc, fontsize=30, color='white')

    norm = ImageNormalize(radio, **norm_params)
    im = ax.imshow(radio, cmap='inferno', origin='lower', norm=norm)

    # colorbar
    #fig.subplots_adjust(top=0.85)
    cbaxes = fig.add_axes([0.127, 0.8, 0.772, 0.02])  # Adjusted colorbar position
    cbar = fig.colorbar(im, cax=cbaxes, ticks=ticks, format='$%.2f$', orientation='horizontal')
    cbar.set_label(r'Surface brightness (mJy beam$^{-1}$)', labelpad=20)
    cbaxes.xaxis.tick_top()
    cbaxes.tick_params(labelsize=20)
    cbaxes.xaxis.label.set_size(25)
    cbaxes.xaxis.set_label_position('top')

    fig.savefig(output_file, bbox_inches='tight')
    plt.close(fig)

# Configuration for L band
l_band_image = '/data/abell_3667/l_band/node35-MFS-I-image-PBcorr.fits'
l_band_norm_params_N = {'vmin': 0.003, 'vmax': 1., 'stretch': SqrtStretch()}
l_band_norm_params_S = {'vmin': 0.005, 'vmax': 1.2, 'stretch': SqrtStretch()}
l_band_ticks_N = [0.03, 0.15, 0.5, 0.8]
l_band_ticks_S = [0.03, 0.3, 1.0]

# Configuration for UHF band
uhf_band_image = '/data/abell_3667/img_uhf/ddcal_005-MFS-image_pbcor.fits'
uhf_band_norm_params_N = {'vmin': 0.003, 'vmax': 4., 'stretch': SqrtStretch()}
uhf_band_norm_params_S = {'vmin': 0.005, 'vmax': 2., 'stretch': SqrtStretch()}
uhf_band_ticks_N = [0.10, 0.8, 2.0, 3.5]
uhf_band_ticks_S = [0.10, 0.6, 1.5]

# Plot configurations for North and South relics
configurations = [
    {
        'image_path': l_band_image,
        'center': [302.77186, -56.4840023], # relic
        #'center': [302.7349515, -55.5070391], # zoom
        'size': [0.57, 0.44], # relic
        #'size': [0.1444700, 0.1550410], # zoom
        'scalebar_kpc': 1000,
        'norm_params': l_band_norm_params_N,
        'ticks': l_band_ticks_N,
        'figsize': (20, 21),
        'output_file': 'radio_relicN_Lband.pdf'
    },
    {
        'image_path': l_band_image,
        'center': [303.5453, -57.05],
        'size': [0.3297, 0.5546],
        'scalebar_kpc': 500,
        'norm_params': l_band_norm_params_S,
        'ticks': l_band_ticks_S,
        'figsize': (10, 22),
        'output_file': 'radio_relicS_Lband.pdf'
    },
    {
        'image_path': uhf_band_image,
        'center': [302.77186, -56.4840023], # relic
        #'center': [302.7349515, -55.5070391], # zoom
        'size': [0.57, 0.44], # relic
        #'size': [0.1444700, 0.1550410], # zoom
        'scalebar_kpc': 1000,
        'norm_params': uhf_band_norm_params_N,
        'ticks': uhf_band_ticks_N,
        'figsize': (20, 21),
        'output_file': 'radio_relicN_UHFband.pdf'
    },
    {
        'image_path': uhf_band_image,
        #'center': [303.5776401, -57.1209859],
        'center': [303.5453, -57.05],
        #'size': [0.57 / 2., 0.6],
        'size': [0.3297, 0.5546],
        'scalebar_kpc': 500,
        'norm_params': uhf_band_norm_params_S,
        'ticks': uhf_band_ticks_S,
        'figsize': (10, 22),
        'output_file': 'radio_relicS_UHFband.pdf'
    }
]

# Generate plots
for config in configurations:
    plot_relic(**config)
