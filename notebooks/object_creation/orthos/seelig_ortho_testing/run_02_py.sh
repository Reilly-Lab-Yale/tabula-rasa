#!/bin/bash
#SBATCH --job-name=fitting_seelig_ortho_cell_type_1
#SBATCH --partition=week
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1    
#SBATCH --mem=16G              
#SBATCH --time=7-00:00:00           
#SBATCH --output=seelig_test_7day_ct1_%j.out    
#SBATCH --error=seelig_test_7day_ct1_%j.err     

module load miniconda
conda activate scmpra 

python 02_ortho_creation_seelig.py
