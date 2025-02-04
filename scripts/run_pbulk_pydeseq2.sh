#!/bin/sh
#SBATCH --job-name=degs   
#SBATCH -A SEYEDAM_LAB           
#SBATCH --nodes=1 
#SBATCH -p highmem
#SBATCH --cpus-per-task=1         
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err
#SBATCH --mem=200G
#SBATCH --array=1  # Array job with 3 tasks (one for each tissue)
#SBATCH --time=1-00:00:00

source ~/miniconda3/bin/activate scanpy_env

# Array of tissues
tissues=('Adrenal')

#'GonadsFemale' 'DiencephalonPituitary' 'Liver' 'CortexHippocampus' 'Heart' 'Gastrocnemius' 'GonadsMale' 'Kidney' 'Adrenal')

# Get the tissue corresponding to the array job task ID
tissue=${tissues[$SLURM_ARRAY_TASK_ID-1]}

# Run the Python script with the current tissue
python3 pseudobulk_pydeseq.py -t "$tissue" -l subtype
