#!/bin/bash
#SBATCH --job-name=fitting_seelig_ortho
#SBATCH --partition=day
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1    
#SBATCH --mem=16G              
#SBATCH --time=06:00:00           
#SBATCH --output=master_%j.out    
#SBATCH --error=master_%j.err     

module load miniconda
conda activate manual 

python 02_ortho_creation_seelig.py
