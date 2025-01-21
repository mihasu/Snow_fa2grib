#!/bin/bash
module load anaconda3
module load gcc
module load openmpi
conda activate epygram

# Directory containing ICMSHHARM files
input_directory="/lustre/tmp/hirlam2/Snowparams_BMM/metcoop_input/"

# Find all ICMSHHARM files in the directory
files=$(find "$input_directory" -type f -name "ICMSHSELE+*.sfx")

# Check if any files were found
if [ -z "$files" ]; then
  echo "No ICMSHSHELE files found in the directory."
  exit 1
fi

# Launch snow_fa2grib.py for each file in parallel using srun
for file in $files; do
  srun --exclusive -n 1 /lustre/tmp/hirlam2/Snowparams_BMM/scripts/snow_fa2grib.py "$file" &
  sleep 2
done

# Wait for all background jobs to finish
wait

echo "All files processed."
