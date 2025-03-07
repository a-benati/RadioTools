#!/usr/bin/env python

import argparse
import os
import subprocess
import sys
import casacore.tables as pt

def run_command(command):
    try:
        result = subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print(result.stdout.decode(), end='')
        print(result.stderr.decode(), end='')
    except subprocess.CalledProcessError as e:
        print(f"Error occurred during command: {command}")
        print(e.stderr.decode())
        sys.exit(1)

def create_point_sources_image(ms_file, scale, size, column, uvmin_l, facet_regions, facet_solutions):
    command = (
        f"wsclean -name {output_dir}/point_sources_uvmin{uvmin_l} "
        f"-j 32 "
        f"-mem 100 "
        f"-abs-mem 100.0 "
        f"-no-update-model-required "
        f"-weight briggs -1 "
        f"-no-mf-weighting "
        f"-size {size} {size} "
        f"-scale {scale} "
        f"-channels-out 15 "
        f"-join-channels "
        f"-nwlayers-factor 3 "
        f"-pol I "
        f"-data-column {column} "
        f"-minuv-l {uvmin_l} "
        f"-niter 1000000 "
        f"-auto-threshold 0.5 "
        f"-auto-mask 2.5 "
        f"-gain 0.1 "
        f"-mgain 0.95 "
        f"-fit-spectral-pol 4 "
        f"-padding 1.3 "
        f"-parallel-deconvolution 1238 "
        f"-multiscale "
        f"-facet-regions {facet_regions} "
        f"-apply-facet-solutions {facet_solutions} amplitude000,phase000 "
        f"{ms_file}"
    )
    run_command(command)

def make_point_sources_mask(uvmin_l):
    command = (
        f"/opt/lofar/DDFacet/SkyModel/./MakeMask.py "
        f"--RestoredIm={output_dir}/point_sources_uvmin{uvmin_l}-MFS-image.fits "
        f"--Th=3.0 "
        f"--Box=100,2"
    )
    run_command(command)

def create_masked_point_sources_image(ms_file, scale, size, column, uvmin_l, facet_regions, facet_solutions):
    command = (
        f"wsclean -name {output_dir}/mask_point_sources_uvmin{uvmin_l} "
        f"-j 32 "
        f"-mem 100 "
        f"-abs-mem 100.0 "
        f"-no-update-model-required "
        f"-weight briggs -1 "
        f"-no-mf-weighting "
        f"-size {size} {size} "
        f"-scale {scale} "
        f"-channels-out 15 "
        f"-join-channels "
        f"-nwlayers-factor 3 "
        f"-pol I "
        f"-data-column {column} "
        f"-minuv-l {uvmin_l} "
        f"-niter 1000000 "
        f"-auto-threshold 0.5 "
        f"-auto-mask 2.5 "
        f"-fits-mask {output_dir}/point_sources_uvmin{uvmin_l}-MFS-image.fits.mask.fits "
        f"-gain 0.1 "
        f"-mgain 0.95 "
        f"-fit-spectral-pol 4 "
        f"-padding 1.3 "
        f"-parallel-deconvolution 1238 "
        f"-multiscale "
        f"-facet-regions {facet_regions} "
        f"-apply-facet-solutions {facet_solutions} amplitude000,phase000 "
        f"{ms_file}"
    )
    run_command(command)

def model_prediction(ms_file, uvmin_l, facet_regions, facet_solutions):
    command = (
        f"wsclean -predict "
        f"-channels-out 15 "
        f"-name {output_dir}/mask_point_sources_uvmin{uvmin_l} "
        f"-facet-regions {facet_regions} "
        f"-apply-facet-solutions {facet_solutions} amplitude000,phase000 "
        f"{ms_file}"
    )
    run_command(command)

def subtract_model_from_uv_plane(ms_file, facet_regions, facet_solutions):
    outcolumn = 'DIFFUSE_SUB'

    ts = pt.table({ms_file}, readonly=False)
    colnames = ts.colnames()

    desc = ts.getcoldesc('DATA')
    desc['name']=outcolumn
    ts.addcols(desc)

    data = ts.getcol('DATA')
    model = ts.getcol('MODEL_DATA')
    ts.putcol(outcolumn,data-model)
    ts.close()

    # command = (
    #     f"wsclean -subtract-model "
    #     f"-channels-out 15 "
    #     f"-name {output_dir}/mask_point_sources_uvmin{uvmin_l} "
    #     f"-facet-regions {facet_regions} "
    #     f"-apply-facet-solutions {facet_solutions} amplitude000,phase000 "
    #     f"{ms_file}"
    # )
    # run_command(command)

def create_extended_sources_image(ms_file, scale, size, uvmin_l, facet_regions, facet_solutions):
    command = (
        f"wsclean -name {output_dir}/sub_uvmin{uvmin_l} "
        f"-j 32 "
        f"-mem 100 "
        f"-abs-mem 100.0 "
        f"-no-update-model-required "
        f"-weight briggs -0.5 "
        f"-no-mf-weighting "
        f"-size {size} {size} "
        f"-scale {scale} "
        f"-channels-out 15 "
        f"-join-channels "
        f"-nwlayers-factor 3 "
        f"-pol I "
        f"-data-column DIFFUSE_SUB "
        f"-niter 1000000 "
        f"-auto-threshold 0.5 "
        f"-auto-mask 2.5 "
        f"-gain 0.1 "
        f"-mgain 0.95 "
        f"-fit-spectral-pol 4 "
        f"-padding 1.3 "
        f"-parallel-deconvolution 1238 "
        f"-multiscale "
        f"-facet-regions {facet_regions} "
        f"-apply-facet-solutions {facet_solutions} amplitude000,phase000 "
        f"{ms_file}"
    )
    run_command(command)

def create_tapered_images(ms_file, scale, size, uvmin_l, facet_regions, facet_solutions, tapers):
    for taper in tapers:
        command = (
            f"wsclean -name {output_dir}/sub_uvmin{uvmin_l}_taper{taper} "
            f"-j 32 "
            f"-mem 100 "
            f"-abs-mem 100.0 "
            f"-no-update-model-required "
            f"-weight briggs -0.5 "
            f"-no-mf-weighting "
            f"-size {size} {size} "
            f"-scale {scale} "
            f"-channels-out 15 "
            f"-join-channels "
            f"-nwlayers-factor 3 "
            f"-pol I "
            f"-data-column DIFFUSE_SUB "
            f"-niter 1000000 "
            f"-auto-threshold 0.5 "
            f"-auto-mask 2.5 "
            f"-gain 0.1 "
            f"-mgain 0.95 "
            f"-fit-spectral-pol 4 "
            f"-padding 1.3 "
            f"-parallel-deconvolution 1238 "
            f"-multiscale "
            f"-facet-regions {facet_regions} "
            f"-apply-facet-solutions {facet_solutions} amplitude000,phase000 "
            f"-taper-gaussian {taper} "
            f"{ms_file}"
        )
        run_command(command)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Point source subtraction in radio images")
    parser.add_argument('--scale', required=True, help='Pixel scale of the image (e.g., 1.9asec)')
    parser.add_argument('--size', required=True, type=int, help='Size of the image')
    parser.add_argument('--column', required=True, help='Column to image (e.g., CORRECTED_DATA)')
    parser.add_argument('--uvmin_l', required=True, type=float, help='UVmin-L parameter')
    parser.add_argument('--facet-regions', required=True, help='Facet region file')
    parser.add_argument('--facet-solutions', required=True, help='Facet solution file (h5 file)')
    parser.add_argument('--tapers', required=True, nargs='+', type=float, help='Tapers to use (e.g., 25 50 75)')
    parser.add_argument('--ms-file', required=True, help='Measurement set file to process')
    args = parser.parse_args()

    # Create output directory
    output_dir = f"/data/abell_3667/sub_sources/temp/sub_sources_uvcut_{args.uvmin_l}"
    os.makedirs(output_dir, exist_ok=True)

    # Execute the functions in sequence
    #create_point_sources_image(args.ms_file, args.scale, args.size, args.column, args.uvmin_l, args.facet_regions, args.facet_solutions)
    #make_point_sources_mask(args.uvmin_l)
    #create_masked_point_sources_image(args.ms_file, args.scale, args.size, args.column, args.uvmin_l, args.facet_regions, args.facet_solutions)
    #model_prediction(args.ms_file, args.uvmin_l, args.facet_regions, args.facet_solutions)
    subtract_model_from_uv_plane(args.ms_file, args.facet_regions, args.facet_solutions)
    create_extended_sources_image(args.ms_file, args.scale, args.size, args.uvmin_l, args.facet_regions, args.facet_solutions)
    create_tapered_images(args.ms_file, args.scale, args.size, args.uvmin_l, args.facet_regions, args.facet_solutions, args.tapers)
