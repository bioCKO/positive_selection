#!/bin/bash
#SBATCH -J rapid_c80TreeBoot
#SBATCH -n 1
#SBATCH -p normal
#SBATCH -o rapid_c80TreeBoot.o%j
#SBATCH -e rapid_c80TreeBoot.e%j
#SBATCH -t 24:00:00
#SBATCH -A tagmap
#SBATCH --mail-type=ALL
#SBATCH --mail-user=grovesdixon@gmail.com

module load launcher





export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
export LAUNCHER_RMI=SLURM
export LAUNCHER_JOB_FILE=rapid_c80TreeBoot

$LAUNCHER_DIR/paramrun
