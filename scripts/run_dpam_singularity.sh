#!/bin/bash
#SBATCH --job-name=dpam
#SBATCH --output=dpam-%j.out
#SBATCH --error=dpam-%j.out
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=128GB

root_dir=$(git rev-parse --show-toplevel)
cd $root_dir

python ./run_dpam_singularity.py \
  --databases_dir databases \
  --input_dir example \
  --dataset test \
  --threads 32 \
  --image_name dpam.sif 