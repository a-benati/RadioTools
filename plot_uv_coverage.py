#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from casacore.tables import table
import argparse

def plot_uv_coverage(ms_file, uvmin=0):
    # Open the MS table
    ms_table = table(ms_file, ack=False)
    uvw = ms_table.getcol('UVW')
    ms_table.close()

    # Compute the uv distances in meters
    u = uvw[:, 0]
    v = uvw[:, 1]
    uv_dist = np.sqrt(u**2 + v**2)

    # Apply uvmin cut
    if uvmin > 0:
        mask = uv_dist >= uvmin
        u = u[mask]
        v = v[mask]

    # Plot the uv coverage
    plt.figure(figsize=(10, 10))
    plt.scatter(u, v, s=1, alpha=0.5, label='UV points', color='b')
    plt.scatter(-u, -v, s=1, alpha=0.5, label='_nolegend_', color='b')  # Include conjugate points for symmetry
    plt.xlabel('u (meters)')
    plt.ylabel('v (meters)')
    plt.title(f'UV Coverage with UVmin = {uvmin} meters')
    plt.legend()
    plt.grid()
    plt.axis('equal')

    # Determine the output file name
    if uvmin > 0:
        output_filename = f'uv_cov_uvmin{uvmin}.pdf'
    else:
        output_filename = 'uv_cov.pdf'

    # Save the plot as a PDF file
    plt.savefig(output_filename, bbox_inches='tight')

    print('UV coverage plot saved as {output_filename}')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot UV coverage from a Measurement Set (MS) file.")
    parser.add_argument("ms_file", help="Path to the Measurement Set (MS) file.")
    parser.add_argument("--uvmin", type=float, default=0, help="Minimum UV distance to include (in meters).")
    args = parser.parse_args()

    plot_uv_coverage(args.ms_file, args.uvmin)
