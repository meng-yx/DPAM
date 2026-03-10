#!/bin/bash
#SBATCH --job-name=precompute_metrics
#SBATCH --array=1-1000%200
#SBATCH --output=logs/precompute_metrics_%A/precompute_metrics-%A_%a.out
#SBATCH --error=logs/precompute_metrics_%A/precompute_metrics-%A_%a.out
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8GB

# -----------------------------------------------------------------------------
# Always run from the repo root
# -----------------------------------------------------------------------------
root_dir=$(git rev-parse --show-toplevel)
cd "$root_dir"


# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
DOMAINOME_DIR="${DOMAINOME_DIR:-./DPAM-AI_AFDB_domainome_v6}"
subset_id=$(printf "%04d" "$SLURM_ARRAY_TASK_ID")

input_csv="${DOMAINOME_DIR}/info/${subset_id}_domains_info.csv"
output_csv="${DOMAINOME_DIR}/metrics/${subset_id}_domains_info.csv"
pdb_dir="${DOMAINOME_DIR}/pdbs"

# Skip if input missing
if [[ ! -f "$input_csv" ]]; then
  echo "Skipping subset ${subset_id}: input not found: $input_csv"
  exit 0
fi

mkdir -p "${DOMAINOME_DIR}/metrics"

echo "Processing subset ${subset_id}"
echo "  input: $input_csv"
echo "  output: $output_csv"
echo "  pdb_dir: $pdb_dir"

# only run if output does not exist
if [[ -f "$output_csv" ]]; then
  echo "Skipping subset ${subset_id}: output already exists: $output_csv"
else
  # Use conda env
  source ~/.bashrc
  conda activate MaSIF

  # Stagger start to avoid all jobs hitting UniProt simultaneously
  sleep $((SLURM_ARRAY_TASK_ID % 60))

  python scripts/python/precompute_metrics.py \
    --input "$input_csv" \
    --output "$output_csv" \
    --pdb-dir "$pdb_dir" \
    --uniprot-timeout 60 \
    --uniprot-max-retries 5
fi
