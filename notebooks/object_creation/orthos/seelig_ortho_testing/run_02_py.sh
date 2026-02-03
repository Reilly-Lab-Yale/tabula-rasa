#!/bin/bash
#SBATCH --job-name=fitting_seelig_ortho
#SBATCH --partition=week
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1    
#SBATCH --mem=16G              
#SBATCH --time=4-00:00:00           
#SBATCH --output=seelig_fit_all_%j.out    
#SBATCH --error=seelig_fit_all_%j.err     

module load miniconda
conda activate scmpra 

python 02_ortho_creation_seelig.py
