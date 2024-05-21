#!/usr/bin/env python3

import sys, os 
from astropy.wcs import WCS
from astropy.visualization import MinMaxInterval, SqrtStretch, ImageNormalize, LogStretch
from astropy.io import fits
import numpy as np
from astropy import units as u
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from lib_fits import flatten
import lib_plot

#image = 'radio/Abell_3667_aFix_pol_I_Farcsec_5pln_cor.fits'
image = 'radio/node34-MFS-image-PBcorr.fits'
header, radio = flatten(image)
wcs=WCS(header)
radio *= 1e3

# north
fig = plt.figure(figsize=(20,21))
ax = fig.add_subplot(1, 1, 1, projection=wcs, slices=('x', 'y'))
lon = ax.coords['ra']
lat = ax.coords['dec']
lon.set_axislabel('Right Ascension (J2000)', size=20)
lat.set_axislabel('Declination (J2000)', size=20)
lon.set_major_formatter('hh:mm:ss')
lat.set_major_formatter('dd:mm')
lat.set_ticklabel(rotation=90) # to turn dec vertical

center = [302.77186,-56.4840023]
size = [0.57/np.cos(center[1]*np.pi/180),0.44]
xlims = [center[0]-size[0]/2., center[0]+size[0]/2.]
ylims = [center[1]-size[1]/2., center[1]+size[1]/2.]
x,y = wcs.wcs_world2pix(xlims*u.deg, ylims*u.deg, 1, ra_dec_order=True)
ax.set_xlim(x[1], x[0])
ax.set_ylim(y[0], y[1])
lon.set_ticklabel(size=20)
lat.set_ticklabel(size=20)

# scalebar
z = 0.0553 # redshift
kpc = 1000 # how many kpc is the scalebar?
lib_plot.add_scalebar(ax, wcs, z, kpc)

norm = ImageNormalize(radio, vmin=0.003, vmax=.6, stretch=SqrtStretch())
#im = ax.imshow(radio, cmap='afmhot', origin='lower', norm=norm)
im = ax.imshow(radio, cmap='inferno', origin='lower', norm=norm)

# colorbar
cbaxes = fig.add_axes([0.127, 0.8, 0.772, 0.02])
cbar = fig.colorbar(im, cax=cbaxes, ticks=[0.03,0.1,0.3,0.5], format='$%.2f$', orientation='horizontal')
cbar.set_label(r'Surface brightness (mJy beam$^{-1}$)', labelpad=20)
cbaxes.xaxis.tick_top()
cbaxes.tick_params(labelsize=20)
cbaxes.xaxis.label.set_size(20)
cbaxes.xaxis.set_label_position('top')

#plt.axis('off')
fig.savefig('radio_relicN.pdf',bbox_inches='tight')

# south
fig = plt.figure(figsize=(10,21))
ax = fig.add_subplot(1, 1, 1, projection=wcs, slices=('x', 'y'))
lon = ax.coords['ra']
lat = ax.coords['dec']
lon.set_axislabel('Right Ascension (J2000)', size=20)
lat.set_axislabel('Declination (J2000)', size=20)
lon.set_major_formatter('hh:mm:ss')
lat.set_major_formatter('dd:mm')
lat.set_ticklabel(rotation=90) # to turn dec vertical

center = [303.5776401,-57.1209859]
size = [0.57/np.cos(center[1]*np.pi/180)/2.,0.44]
xlims = [center[0]-size[0]/2., center[0]+size[0]/2.]
ylims = [center[1]-size[1]/2., center[1]+size[1]/2.]
x,y = wcs.wcs_world2pix(xlims*u.deg, ylims*u.deg, 1, ra_dec_order=True)
ax.set_xlim(x[1], x[0])
ax.set_ylim(y[0], y[1])
lon.set_ticklabel(size=20)
lat.set_ticklabel(size=20)

# scalebar
z = 0.0553 # redshift
kpc = 500 # how many kpc is the scalebar?
lib_plot.add_scalebar(ax, wcs, z, kpc)

norm = ImageNormalize(radio, vmin=0.005, vmax=.4, stretch=SqrtStretch())
im = ax.imshow(radio, cmap='inferno', origin='lower', norm=norm)

# colorbar
cbaxes = fig.add_axes([0.127, 0.8, 0.772, 0.02])
cbar = fig.colorbar(im, cax=cbaxes, ticks=[0.03,0.1,0.3,0.5], format='$%.2f$', orientation='horizontal')
cbar.set_label(r'Surface brightness (mJy beam$^{-1}$)', labelpad=20)
cbaxes.xaxis.tick_top()
cbaxes.tick_params(labelsize=20)
cbaxes.xaxis.label.set_size(25)
cbaxes.xaxis.set_label_position('top')

fig.savefig('radio_relicS.pdf',bbox_inches='tight')
