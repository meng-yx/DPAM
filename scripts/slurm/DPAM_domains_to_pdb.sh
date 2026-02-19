#!/bin/bash
#SBATCH --job-name=DPAM_domains_to_pdb
#SBATCH --array=961
#SBATCH --output=logs/DPAM_domains_to_pdb_%A/domains_to_pdb-%A_%a.out
#SBATCH --error=logs/DPAM_domains_to_pdb_%A/domains_to_pdb-%A_%a.out
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=12GB

# -----------------------------------------------------------------------------
# Absolute paths - configure according to your environment
# -----------------------------------------------------------------------------
CHUNKS_BASE="${CHUNKS_BASE:-/scratch/diazrovi/DPAM/processing_human_afdb_v6/chunks_human_afdb}"
DOMAINOME_DIR="${DOMAINOME_DIR:-./DPAM-AI_AFDB_domainome_v6}"


# -----------------------------------------------------------------------------
# Always run from the repo root
root_dir=$(git rev-parse --show-toplevel)
cd "$root_dir"

# Use conda env
source ~/.bashrc
conda activate MaSIF


# Zero-padded subset_id (0001 to 1000)
subset_id=$(printf "%04d" "$SLURM_ARRAY_TASK_ID")

domains="${CHUNKS_BASE}/job_${subset_id}/${subset_id}_domains"
in_pdb_dir="${CHUNKS_BASE}/job_${subset_id}/${subset_id}"
out_pdb_dir="${DOMAINOME_DIR}/pdbs"
out_csv="${DOMAINOME_DIR}/info/${subset_id}_domains_info.csv"

# Skip if domains file does not exist
if [[ ! -f "$domains" ]]; then
  echo "Skipping subset ${subset_id}: domains file not found: $domains"
  exit 0
fi

# Create output directories
mkdir -p "$out_pdb_dir"
mkdir -p "$(dirname "$out_csv")"

echo "Processing subset ${subset_id}"
echo "  domains: $domains"
echo "  in_pdb_dir: $in_pdb_dir"
echo "  out_pdb_dir: $out_pdb_dir"
echo "  out_csv: $out_csv"

# Run the script
python ./scripts/python/DPAM_domains_to_pdb.py \
  --domains "$domains" \
  --in_pdb_dir "$in_pdb_dir" \
  --out_pdb_dir "$out_pdb_dir" \
  --out_csv "$out_csv"
