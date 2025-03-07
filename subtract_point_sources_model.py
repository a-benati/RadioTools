#!/usr/bin/env python3

import casacore.tables as pt

ms = '/data/abell_3667/msfiles/1685906777_sdp_l0-Abell3667-corr_avg_DP3_copy.ms'

ts = pt.table(ms, readonly=False)

data = ts.getcol('DATA')
model = ts.getcol('MODEL_DATA')

subtracted_data = data - model
# Print sizes and shapes
# print(f"DATA column size: {data.shape}")
# print(f"MODEL_DATA column size: {model.shape}")
ts.putcol('DIFFUSE_SUB', subtracted_data)
ts.close()