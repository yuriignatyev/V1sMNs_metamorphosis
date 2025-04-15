#!/bin/bash
#
#----------------------------------
# CellBender on GPU cluster
#----------------------------------

#SBATCH --job-name=cellbender_gpu_example
#SBATCH --output=cellbender_output_log.txt
#SBATCH --ntasks=20
#SBATCH --time=36:00:00
#SBATCH --mem=99G
#SBATCH --mail-user=your_email@institute.edu
#SBATCH --mail-type=ALL
#SBATCH --no-requeue
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --constraint=GTX1080Ti
#SBATCH --export=NONE

unset SLURM_EXPORT_ENV
export OMP_NUM_THREADS=20

# Load CUDA and conda environment with CellBender installed
module load cuda
module load miniforge3/24.7.1
conda activate cellbender2

# set paths for raw input from CellRanger count and your output
RAW_INPUT="/path/to/your/sample/raw_feature_bc_matrix.h5"
OUTPUT="/path/to/output/output_sample"

# run CellBender (here multiple fpr's are used, adjust as needed. In the end, outputs from fpr = 0.01 were used for further downstream)
cellbender remove-background \
    --input "$RAW_INPUT" \
    --output "$OUTPUT" \
    --epochs 150 \
    --cuda \
    --fpr 0.01 0.05 0.1 0.3 \
    --empty-drop-training-fraction 0.3

conda deactivate
