#!/bin/bash
#SBATCH --job-name=phastest_master
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=7-00:00:00

#SBATCH --output=xtt_project/mobilome/prophages/phastest/logs/phastest_master_%j.out
#SBATCH --error=xtt_project/mobilome/prophages/phastest/logs/phastest_master_%j.err

echo "================================="
echo "PHASTEST MASTER PIPELINE"
echo "Start: $(date)"
echo "Node: $(hostname)"
echo "================================="

cd ~/xtt_population_genomics

source ~/xtt_population_genomics/envs/xtt_mobilome_clean/bin/activate

which python

python -u \
xtt_project/mobilome/prophages/phastest/scripts/phastest_master_pipeline.py

echo "================================="
echo "Finished: $(date)"
echo "================================="