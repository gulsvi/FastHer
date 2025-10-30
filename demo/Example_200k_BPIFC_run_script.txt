#!/bin/bash
#
#SBATCH --job-name=Example_BPIFC
#SBATCH --output=Example_200k_BPIFC_run.log

srun R --vanilla < Example_200k_BPIFC_run.R

