#!/bin/bash

# Set the environment variables for an inter-annual run
# just putting in the variables which we may want to make configurable
# we can always add more here if we want

# source code
export WW3_EXE_DIR=/home/$USER/code/WW3/model/exe/

# directories corresponding to this configuration
CONFIG_DIR=/home/gfearon/code/somisana-ww3/configs/sa_west_02
export RUN_DIR=$CONFIG_DIR/RUN_01
export GRID_DIR=$CONFIG_DIR/GRID/
export BRY_DIR=$CONFIG_DIR/SPEC_CMEMS/
export SRF_DIR=/mnt/c/Users/GilesF/Downloads/ERA5/

# MPI settings
export MPI_NUM_PROCS=4

# time period for interannual run
MONTH_START="2009-01"
MONTH_END="2009-02"
