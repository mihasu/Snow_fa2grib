#!/bin/bash
module load anaconda3
module load gcc
module load openmpi
conda activate epygram

# Directory containing ICMSHHARM files
input_directory="/lustre/tmp/hirlam2/Snowparams_BMM/metcoop_input/"

# Find all ICMSHHARM files in the directory
files=$(find "$input_directory" -type f -name "ICMSHHARM+*.sfx")

# Check if any files were found
if [ -z "$files" ]; then
  echo "No ICMSHHARM files found in the directory."
  exit 1
fi

# Launch snow_fa2grib.py for each file in parallel using srun
for file in $files; do
  srun --exclusive -n 1 /lustre/tmp/hirlam2/Snowparams_BMM/scripts/Snow_fa2grib/snow_fa2grib.py "$file" &
done

# Wait for all background jobs to finish
wait

echo "All files processed."
