#!/bin/bash
#SBATCH --job-name=download_database
#SBATCH --output=slurm-%A.out  
#SBATCH --error=slurm-%A.out 
#SBATCH --time=12:00:00  
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G

root_dir=$(git rev-parse --show-toplevel)
cd $root_dir

echo "Downloading databases..."
wget -c --no-check-certificate https://conglab.swmed.edu/DPAM/databases.tar.gz

echo "Decompressing databases..."
pigz -p $SLURM_CPUS_PER_TASK -dc databases.tar.gz | tar xf -

echo "Pulling singularity image..."
singularity pull dpam.sif docker://conglab/dpam

echo "Done"