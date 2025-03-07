#!/usr/bin/env python3

import numpy as np

c = 3e8

nu_max_L = 1682.1156e6 # max L band
nu_max_UHF = 1088e6 # max UHF band
nu_min_L = 886.303e6 # min L band
nu_min_UHF = 544e6 # min UHF band
nu_c_L = 1280e6 # central L band
nu_c_UHF = 816e6 # central UHF band
nu_c_both = 1113e6 # central both bands

band = "UHF"  # Può essere "UHF", "L" o "both"

if band == "UHF":
    nu_max = nu_max_UHF
    nu_min = nu_min_UHF
    nu_c = nu_c_UHF
elif band == "L":
    nu_max = nu_max_L
    nu_min = nu_min_L
    nu_c = nu_c_L
elif band == "both":
    nu_max = nu_max_L
    nu_min = nu_min_UHF
    nu_c = nu_c_both
else:
    raise ValueError("Banda non valida! Usa 'UHF', 'L' o 'both'.")

delta_nu = nu_max - nu_min

nu_max_2 = nu_max**2
nu_min_2 = nu_min**2
delta_nu_2 = nu_max_2 - nu_min_2

l_max = c / nu_min
l_min = c / nu_max
delta_l = l_max - l_min
l_max_2 = l_max**2
l_min_2 = l_min**2
delta_l_2 = l_max_2 - l_min_2

dnu_chan = 8.5e6
dl_chan_2 = (2*c**2*dnu_chan)/(nu_c**3)*(1+0.5*((dnu_chan)/(nu_c))**2)

dphi = 2.*np.sqrt(3.) / delta_l_2
max_scale = np.pi / l_min_2
phi_max = np.sqrt(3) / dl_chan_2

print("dphi = ", dphi)
print("max scale = ", max_scale)
print("phi_max = ", phi_max)