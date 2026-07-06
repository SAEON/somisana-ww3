#!/bin/bash

# This script is designed to run on the chpc
# it submits a series of run_inter.pbs job scripts, each of which is for a few months
# this is needed to ensure we stay under the 48 hr job length for each job submission

CURRENT_DATE="2008-12-01"
END_DATE="2013-12-31"
RSTFLAG=0 # set to 1 if you want to restart from a previous months output
prev_jobid=""

# we need to give the absolute path to the script, even though we are running this bash script from the same directory!
SCRIPT="/home/gfearon/lustre/somisana-ww3/configs/sa_west_02/RUN_01/run_inter.pbs"

while [ "$(date -d "$CURRENT_DATE" +%Y%m)" -le "$(date -d "$END_DATE" +%Y%m)" ]; do

    # using jobs of 4 month increments, so add 3 months to get the end dates for this job
    END_DATE_JOB=$(date -d "$CURRENT_DATE +3 month" +%Y-%m-01)

    MONTH_START=$(date -d "$CURRENT_DATE" +%Y-%m)
    MONTH_END=$(date -d "$END_DATE_JOB" +%Y-%m)

    if [ -z "$prev_jobid" ]; then
        # First job (no dependency)
        jobid=$(qsub -v MONTH_START=$MONTH_START,MONTH_END=$MONTH_END,RSTFLAG=$RSTFLAG $SCRIPT)
    else
        # include dependency of previous job
        jobid=$(qsub -W depend=afterok:$prev_jobid -v MONTH_START=$MONTH_START,MONTH_END=$MONTH_END,RSTFLAG=$RSTFLAG $SCRIPT)
    fi

    echo "Submitted run from $MONTH_START to $MONTH_END as job $jobid"
    echo
    prev_jobid=$(echo $jobid | cut -d. -f1)  # Extract numeric job ID

    # Advance by one month
    CURRENT_DATE=$(date -d "$END_DATE_JOB +1 month" +%Y-%m-01)

    # Always set restart to 1 after the first job
    RSTFLAG=1

done
