#!/bin/bash
#SBATCH --job-name=dpam
#SBATCH --output=logs/dpam_array_%A/dpam-%A_%a.out
#SBATCH --error=logs/dpam_array_%A/dpam-%A_%a.out
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=128GB

# Always run from the root directory
root_dir=$(git rev-parse --show-toplevel)
cd $root_dir

# Print the resource usage
echo "Resource usage:"
echo "Memory: $SLURM_MEM_PER_NODE"
echo "CPUs: $SLURM_CPUS_PER_TASK"


# Define path to the _struc.list file containing the UniProt IDs of the models to download
list_path=pacini_results/subset_${SLURM_ARRAY_TASK_ID}_struc.list
list_path=$(realpath $list_path)
list_dir=$(dirname $list_path)
list_name=$(basename "$list_path" "_struc.list")
echo "List name: $list_name"

# Prepare the input files
python ./scripts/python/prepare_input.py \
  --list_path $list_path 

# Start running DPAM-AI
start_time=$(date '+%Y-%m-%d %H:%M:%S')
echo "$start_time - Executing run_dpam_singularity.py..."
python ./run_dpam_singularity.py \
  --databases_dir databases \
  --input_dir $list_dir \
  --dataset $list_name \
  --threads $SLURM_CPUS_PER_TASK \
  --image_name dpam.sif 

end_time=$(date '+%Y-%m-%d %H:%M:%S')
echo "$end_time - Done"

# Calculate the time taken
start_sec=$(date -d "$start_time" +%s)
end_sec=$(date -d "$end_time" +%s)
duration=$((end_sec - start_sec))
printf "Time taken: %02d:%02d:%02d\n" $((duration/3600)) $(( (duration%3600)/60 )) $((duration%60))