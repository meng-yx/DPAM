#!/bin/bash
#SBATCH --job-name=jupyter
#SBATCH --partition=standard
#SBATCH --cpus-per-task=64
#SBATCH --mem=384G
#SBATCH --time=02:00:00
#SBATCH --output=logs/jupyter-%j.out
#SBATCH --error=logs/jupyter-%j.out



echo "Starting Jupyter job on $(hostname)"
echo "SLURM_JOB_ID=$SLURM_JOB_ID"


# --- pick a port ---
PORT=$(python -c "
import socket
s=socket.socket()
s.bind(('',0))
print(s.getsockname()[1])
s.close()
")

# --- activate conda and launch Jupyter ---
set +u
source ~/miniconda3/etc/profile.d/conda.sh
conda activate MaSIF
set -u

jupyter notebook \
  --no-browser \
  --ip=0.0.0.0 \
  --port=${PORT} \
  --ServerApp.allow_remote_access=True 