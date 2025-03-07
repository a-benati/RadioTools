#!/usr/bin/env python3

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

# Apri i file FITS
q_hdul = fits.open("/data/abell_3667/calibration/polcal/A3667/pyrmsinth/Lband/output/_di_q.fits")
u_hdul = fits.open("/data/abell_3667/calibration/polcal/A3667/pyrmsinth/Lband/output/_di_u.fits")

q_map = q_hdul[0].data  # Parte reale (Q)
u_map = u_hdul[0].data  # Parte immaginaria (U)
hdr = q_hdul[0].header  # Header per gli assi
q_hdul.close()
u_hdul.close()

# Coordinate del pixel di interesse (esempio: centro)
x1, y1 = 1404, 1700
x2, y2 = 1405, 1700
q_spectrum1 = q_map[:, y1, x1]  
u_spectrum1 = u_map[:, y1, x1]  

# Calcola P(ϕ) = sqrt(Q² + U²)
p_spectrum1 = np.sqrt(q_spectrum1**2 + u_spectrum1**2)

q_spectrum2 = q_map[:, y2, x2]  
u_spectrum2 = u_map[:, y2, x2]  

# Calcola P(ϕ) = sqrt(Q² + U²)
p_spectrum2 = np.sqrt(q_spectrum2**2 + u_spectrum2**2)

# Estrarre asse della Faraday Depth (ϕ)
n_phi = hdr["NAXIS3"]
phi_min = hdr["CRVAL3"]
delta_phi = hdr["CDELT3"]
phi_axis = phi_min + delta_phi * np.arange(n_phi)

# Trova tutti i picchi
peaks1, properties1 = find_peaks(np.abs(p_spectrum1), height=0)
peaks2, properties2 = find_peaks(np.abs(p_spectrum2), height=0)

# Ordina i picchi per altezza e prendi i due maggiori
sorted_indices1 = np.argsort(properties1["peak_heights"])[::-1]  # Ordina in ordine decrescente
top_two_peaks1 = peaks1[sorted_indices1[:2]]
sorted_indices2 = np.argsort(properties2["peak_heights"])[::-1]  # Ordina in ordine decrescente
top_two_peaks2 = peaks2[sorted_indices2[:2]]

# Stampiamo le profondità di Faraday dei due picchi maggiori
print(f"Due picchi principali pixel 1 a: {phi_axis[top_two_peaks1]} rad/m²")
print(f"Due picchi principali pixel 2 a: {phi_axis[top_two_peaks2]} rad/m²")

# peak_indices = np.argsort(np.abs(p_spectrum1))[-2:]  # Trova i due picchi principali  
# peak_phi_values = phi_axis[peak_indices]  
# print(f"Faraday Depth peaks: {peak_phi_values}")

# Plotta lo spettro di Faraday
# plt.figure(figsize=(8, 5))
# plt.plot(phi_axis, p_spectrum, label="|FDF(ϕ)|", color="black")
# plt.plot(phi_axis, q_spectrum, label="Re(FDF) (Q)", linestyle="dashed")
# plt.plot(phi_axis, u_spectrum, label="Im(FDF) (U)", linestyle="dotted")
# plt.xlabel("Faraday Depth ϕ [rad/m²]")
# plt.ylabel("Polarized Intensity")
# #plt.xlim(0, 100)
# plt.legend()
# plt.show()
plt.figure(figsize=(8,5))
plt.plot(phi_axis, np.abs(p_spectrum1), label="|FDF(ϕ)| - Pixel 1")
plt.plot(phi_axis, np.abs(p_spectrum2), label="|FDF(ϕ)| - Pixel 2", linestyle="dashed")
plt.xlabel("Faraday Depth ϕ [rad/m²]")
plt.ylabel("|FDF(ϕ)|")
plt.legend()
plt.show()
plt.savefig("/data/abell_3667/calibration/polcal/A3667/pyrmsinth/Lband/output/faraday_spectrum_pixels.png")
