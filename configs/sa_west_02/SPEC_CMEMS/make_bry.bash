#!/bin/bash

# with the wavespectra environment activated!

python /home/gfearon/code/somisana-ww3/cli.py make_bry_cmems_monthy \
	--input_dir /mnt/c/Users/GilesF/Downloads/CMEMS_WAV/eez/ \
	--lon_file /home/gfearon/code/somisana-ww3/configs/sa_west_02/GRID/lon.dat \
	--lat_file /home/gfearon/code/somisana-ww3/configs/sa_west_02/GRID/lat.dat \
	--mask_file /home/gfearon/code/somisana-ww3/configs/sa_west_02/GRID/mask.dat \
	--output_dir /home/gfearon/code/somisana-ww3/configs/sa_west_02/SPEC_CMEMS \
	--month_start 2009-01 \
	--month_end 2009-02 

