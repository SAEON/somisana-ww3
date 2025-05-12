#!/bin/bash

#unalias cp
#unalias mv
#limit coredumpsize unlimited
CP=/bin/cp
MV=/bin/mv
LN=/bin/ln
MK=/bin/mkdir

# read in the environment variables
source myenv_inter.sh

# create the scratch and output dirs
# scratch is where we'll run the model and output is where we'll send the files we want to keep
$MK $RUN_DIR/scratch
$MK $RUN_DIR/output
cd $RUN_DIR/scratch

# grid files
$CP $WW3_EXE_DIR/ww3_grid .
$CP $GRID_DIR/*.dat .
$CP $RUN_DIR/ww3_grid.nml .
$CP $RUN_DIR/namelists_config.nml .
#./ww3_grid | tee ww3_grid.out

# boundary files
$CP $WW3_EXE_DIR/ww3_bounc .
$CP $RUN_DIR/ww3_bounc.nml .

# surface files
$CP $WW3_EXE_DIR/ww3_prnc .
$CP $RUN_DIR/ww3_prnc.nml .

# shel
$CP $WW3_EXE_DIR/ww3_shel .
$CP $RUN_DIR/ww3_shel.nml ww3_shel_template.nml

# postprocessing script to create nc outputs
$CP $WW3_EXE_DIR/ww3_ounf .
$CP $RUN_DIR/ww3_ounf.nml ww3_ounf_template.nml

CURRENT_DATE="${MONTH_START}-01"
END_DATE="${MONTH_END}-01"

while [ "$(date -d "$CURRENT_DATE" +%Y%m)" -le "$(date -d "$END_DATE" +%Y%m)" ]; do
    
    # Format months
    MONTH_NOW=$(date -d "$CURRENT_DATE" +%Y_%m) 
    MONTH_NEXT=$(date -d "$CURRENT_DATE +1 month" +%Y_%m)
    
    echo "MONTH_NOW:  $MONTH_NOW"
    echo "MONTH_NEXT: $MONTH_NEXT"
    echo

    # create the spectral boundary conditions
    $CP $BRY_DIR/$MONTH_NOW/spec.list .
    #./ww3_bounc | tee ww3_bounc_$MONTH_NOW.out

    # create the surface forcing
    # we want to create a single file called wind.nc which includes u10 and v10 and extends to the first day of the following month so that we have full coverage over
    MONTH_NOW_ERA5=$(date -d "$CURRENT_DATE" +Y%YM%m) 
    MONTH_NEXT_ERA5=$(date -d "$CURRENT_DATE +1 month" +Y%YM%m)
    ufile=$SRF_DIR/ERA5_ecmwf_U10_$MONTH_NOW_ERA5.nc
    vfile=$SRF_DIR/ERA5_ecmwf_V10_$MONTH_NOW_ERA5.nc
    ufile_next=$SRF_DIR/ERA5_ecmwf_U10_$MONTH_NEXT_ERA5.nc
    vfile_next=$SRF_DIR/ERA5_ecmwf_V10_$MONTH_NEXT_ERA5.nc
    # 
    # copy over the wind files and make time a record (unlimited) dimension
    ncks --mk_rec_dmn time $ufile wind_u.nc
    ncks --mk_rec_dmn time $vfile wind_v.nc
    # extract the first time-step of the next month and make time a record (unlimited) dimension
    ncks --mk_rec_dmn time -d time,0,0 $ufile_next wind_u_next.nc
    ncks --mk_rec_dmn time -d time,0,0 $vfile_next wind_v_next.nc
    # and append the first time step of next month to the winds for this month
    ncrcat wind_u.nc wind_u_next.nc wind_u_extended.nc
    ncrcat wind_v.nc wind_v_next.nc wind_v_extended.nc
    # now combine the two wind components into a single file
    cp wind_u_extended.nc wind.nc
    ncks -A wind_v_extended.nc wind.nc
    # and clean up
    rm wind_u* wind_v*

    # run the model
    YMD_FIRSTDAY=$(date -d "$CURRENT_DATE" +%Y%m%d)
    YMD_NEXT_FIRSTDAY=$(date -d "$CURRENT_DATE +1 month" +%Y%m%d)
    YMD_LASTDAY=$(date -d "$CURRENT_DATE +1 month -1 day" +%Y%m%d)

    sed -e 's/YMD_FIRSTDAY/'$YMD_FIRSTDAY'/g' -e 's/YMD_NEXT_FIRSTDAY/'$YMD_NEXT_FIRSTDAY'/g' -e 's/YMD_LASTDAY/'$YMD_LASTDAY'/g' < ww3_shel_template.nml > ww3_shel.nml
    #mpirun -np $MPI_NUM_PROCS ./ww3_shel >& ww3_shel_$MONTH_NOW.out

    # create netcdf from the ww3 binary files
    sed -e 's/YMD_FIRSTDAY/'$YMD_FIRSTDAY'/g' < ww3_ounf_template.nml > ww3_ounf.nml
    #./ww3_ounf | tee ww3_ounf_$MONTH_NOW.out
    
    # move the output
    $MV ww3.$(date -d "$CURRENT_DATE" +%Y%m).nc ../output
    $CP restart.ww3 ../output/restart.$(date -d "$CURRENT_DATE" +%Y%m).ww3

    # Advance by one month
    CURRENT_DATE=$(date -d "$CURRENT_DATE +1 month" +%Y-%m-01)
done

