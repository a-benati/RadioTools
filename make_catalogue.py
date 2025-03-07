#!/usr/bin/env python3

import os, sys, argparse, logging
import numpy as np
from astropy.io import fits as pyfits
from astropy.wcs import WCS as pywcs
from astropy.coordinates import match_coordinates_sky
from astropy.coordinates import SkyCoord
import astropy.units as u
import pyregion
import bdsf
import lsmtool as lsm
from astropy.table import Table
sys.path.append('/homes/a.benati/packages/scripts/')
from lib_fits import AllImages

parser = argparse.ArgumentParser(description='Make spectral index maps, e.g. spidxmap.py --region ds9.reg --noise --sigma 5 --save *fits')
parser.add_argument('images', nargs='+', help='List of images to use for spidx')

args = parser.parse_args()

all_images = AllImages(args.images)

img_cat = all_images[0].imagefile+'.cat'
if not os.path.exists(img_cat):
    bdsf_img = bdsf.process_image(all_images[0].imagefile, rms_box=(100,10), \
                                          thresh_pix=5, thresh_isl=3, atrous_do=False, \
                                          adaptive_rms_box=True, adaptive_thresh=100, rms_box_bright=(30,10), quiet=True)
    # img = bdsf.process_image(image_name, rms_box=rmsbox, \
    #     thresh_isl=int(threshisl), atrous_do=atrous_do, atrous_jmax=3, \
    #     adaptive_rms_box=True, adaptive_thresh=100, rms_box_bright=(30,10), stop_at=stop_at, quiet=True, debug=False)


    bdsf_img.write_catalog(outfile=img_cat, catalog_type='srl', format='fits', clobber=True)
    bdsf_img.write_catalog(outfile=img_cat.replace('.cat', '.skymodel'), catalog_type='gaul', format='bbs', bbs_patches='source', clobber=True, srcroot='src')
else:
    logging.warning('%s already exists, using it.' % img_cat)

cat = Table.read(img_cat)
# remove extended sources
extended_src = (cat['Peak_flux'] / cat['Total_flux']) < 0.1 # ~extended source
extended_src[cat['S_Code'] == 'M'] = True # multiple-gaussian source
extended_src[cat['S_Code'] == 'C'] = True # one gaussian + other sources island
# remove same sources from skymodel
cat_lsm = lsm.load(img_cat.replace('.cat', '.skymodel'))
# TODO this spams the logging
for srcid in cat[extended_src]['Source_id']:
    cat_lsm.remove(f'Patch == src_patch_s{srcid}')
cat.remove_rows(np.argwhere(extended_src))
all_images[0].cat = cat
all_images[0].cat_lsm = cat_lsm
logging.debug('%s: Number of sources detected: %i; removed %i extended sources.' % (all_images[0].imagefile, len(all_images[0].cat), sum(extended_src)) )
