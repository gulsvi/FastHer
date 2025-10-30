#!/bin/bash
#
#SBATCH --job-name=Example_HIRA
#SBATCH --output=Example_200k_HIRA_run.log

srun R --vanilla < Example_200k_HIRA_run.R
