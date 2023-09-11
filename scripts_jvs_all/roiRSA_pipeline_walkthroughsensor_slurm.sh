#!/bin/bash
#SBATCH -A ACLARKE-SL3-CPU
#SBATCH -p cclake
#SBATCH -N 1
#SBATCH -n 12
#SBATCH -t 12:00:00
#SBATCH -J matlab-test
#SBATCH -e matlab_test.err

module purge
module load matlab

matlab -nodisplay -r "addpath('/home/ac584/rds/rds-semconnect-6yHdsDfiMLk/software/spm12/'); spm eeg; addpath('/home/ac584/rds/rds-semconnect-6yHdsDfiMLk/MEG_objects/OscRSA_2016/scripts/scripts_jvs_all/'); roiRSA_pipeline_walkthroughsensor; quit"
