#!/usr/bin/env python3

import argparse
import os
import subprocess
import sys
import casacore.tables as pt
import numpy as np

def get_resolution(ms):
   uvmax = get_uvwmax(ms)
   t = pt.table(ms + '/SPECTRAL_WINDOW', ack=False)
   freq = np.median(t.getcol('CHAN_FREQ'))
   print('Central freq [MHz] = ', freq/1e6)
   print('Longest baselines [km] = ', uvmax/1e3)
   t.close()
   res = 1.22*3600.*180.*((299792458./freq )/uvmax)/np.pi
   print('Resolution [arcsec] = ', res)
   return None

  
def get_uvwmax(ms):
    ''' Find the maximum squared sum of UVW coordinates.
    
    Args:
        ms (str): path to a Measurement Set.
    Returns:
        None
    '''
    t = pt.table(ms, ack=False)
    uvw = t.getcol('UVW')
    ssq = np.sqrt(np.sum(uvw**2, axis=1))
    print(uvw.shape)
    t.close()
    return np.max(ssq)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Get some statistics from a ms file.")
    parser.add_argument('--res', dest='res', action='store_true', help='Print the resolution of the input ms files')
    parser.add_argument('ms', nargs='+', help='List of input ms files.')

    args = parser.parse_args()
    if args.res:
        get_resolution(*args.ms)
