#!/bin/bash

# Set the environment variables for an inter-annual run
# just putting in the variables which we may want to make configurable
# we can always add more here if we want

# source code
export WW3_EXE_DIR=/home/$USER/WW3/model/exe/

# directories corresponding to this configuration
CONFIG_DIR=/home/gfearon/lustre/somisana-ww3/configs/sa_west_02
export RUN_DIR=$CONFIG_DIR/RUN_01
export GRID_DIR=$CONFIG_DIR/GRID/
export BRY_DIR=$CONFIG_DIR/SPEC_CMEMS/
export SRF_DIR=/home/gfearon/lustre/DATA/ERA5/eez

# MPI settings
export MPI_NUM_PROCS=120

# time period for interannual run
MONTH_START="2010-01"
MONTH_END="2010-12"
