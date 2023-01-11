#!/bin/bash

export PBS_JOBID="example_job"
export PBS_ARRAY_INDEX=1
export N_CORES=1
export N_RANGER_THREADS=1
export VARIMP_RUN=FALSE
export N_PERMS=10 # max 1,000

R -f ../scripts/run_rf_experiment.R