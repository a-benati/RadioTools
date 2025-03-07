#!/usr/bin/env python

import casacore.tables as pt
import os

def copy_ms_without_model_data(original_ms, new_ms):
    # Open the original table
    with pt.table(original_ms, readonly=True) as ts:
        colnames = ts.colnames()

        # If 'MODEL_DATA' exists, copy all columns except 'MODEL_DATA'
        if 'MODEL_DATA' in colnames:
            print("'MODEL_DATA' column exists. Creating a new MS without 'MODEL_DATA'...")
            
            # Create the new table description excluding 'MODEL_DATA'
            new_desc = {}
            for colname in colnames:
                if colname != 'MODEL_DATA':
                    new_desc[colname] = ts.getcoldesc(colname)
            
            # Get the data manager information
            dm_info = ts.getdminfo()

            # Remove the 'MODEL_DATA' entry from the data manager info if it exists
            if 'MODEL_DATA' in dm_info:
                del dm_info['MODEL_DATA']
            
            # Create the new table without 'MODEL_DATA' and copy the Data Manager info
            new_table = pt.table(new_ms, new_desc, dm_info=dm_info, nrow=ts.nrows(), readonly=False)

            # Copy the valid columns to the new MS
            for colname in colnames:
                if colname != 'MODEL_DATA':
                    print(f"Copying column: {colname}")
                    new_table.putcol(colname, ts.getcol(colname))

            new_table.close()
            print("New MS created without 'MODEL_DATA'.")
        else:
            print("'MODEL_DATA' column does not exist. No need to remove.")

# Path to the original MS and the new MS without the corrupted 'MODEL_DATA'
original_ms = '/data/abell_3667/msfiles/1685906777_sdp_l0-Abell3667-corr_avg_DP3_copy.ms'
new_ms = '/data/abell_3667/msfiles/1685906777_sdp_l0-Abell3667-corr_avg_DP3_no_model_data.ms'

# Check if the new MS already exists and remove it if necessary
if os.path.exists(new_ms):
    os.remove(new_ms)

# Create the new MS without 'MODEL_DATA'
copy_ms_without_model_data(original_ms, new_ms)
