#!/bin/bash
#SBATCH --job-name=dpam
#SBATCH --output=logs/dpam-%j.out
#SBATCH --error=logs/dpam-%j.out
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
list_path="$1"
list_path=$(realpath $list_path)
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
  --input_dir input \
  --dataset $list_name \
  --threads 32 \
  --image_name dpam.sif 

end_time=$(date '+%Y-%m-%d %H:%M:%S')
echo "$end_time - Done"

# Calculate the time taken
start_sec=$(date -d "$start_time" +%s)
end_sec=$(date -d "$end_time" +%s)
duration=$((end_sec - start_sec))
printf "Time taken: %02d:%02d:%02d\n" $((duration/3600)) $(( (duration%3600)/60 )) $((duration%60))