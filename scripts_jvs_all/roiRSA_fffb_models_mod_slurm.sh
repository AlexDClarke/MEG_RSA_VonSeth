#!/bin/bash
#SBATCH -A ACLARKE-SL3-CPU
#SBATCH -p cclake
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --time=12:00:00
#SBATCH -J matlab-test
#SBATCH -e matlab_test.err

module purge
module load matlab

matlab -nodisplay -r "addpath('/rds/project/rds-6yHdsDfiMLk/MEG_objects/OscRSA_2016/scripts/scripts_jvs_all'); roiRSA_fffb_mod_par; quit"

