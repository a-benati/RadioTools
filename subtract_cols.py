#!/usr/bin/env python3

from casacore.tables import table

# MS files
ms_file_original = '/data/abell_3667/msfiles/1685906777_sdp_l0-Abell3667-corr_avg_DP3.ms'
ms_file_copy = '/data/abell_3667/msfiles/1685906777_sdp_l0-Abell3667-corr_avg_DP3_copy.ms'

# Open the two MS files
with table(ms_file_original, readonly=False) as original_ms:
    original_model_data = original_ms.getcol('MODEL_DATA')

with table(ms_file_copy, readonly=True) as copy_ms:
    mask_model_data = copy_ms.getcol('MODEL_DATA')

# Subtraction between the two MODEL_DATA
new_model_data = original_model_data - mask_model_data

# Save the result in the MODEL_DATA of the original MS
with table(ms_file_original, readonly=False) as original_ms:
    original_ms.putcol('MODEL_DATA', new_model_data)

print("Subtraction complete.")
