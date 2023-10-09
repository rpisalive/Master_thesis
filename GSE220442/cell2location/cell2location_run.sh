#!/bin/bash

#SBATCH --ntasks-per-node=8
#SBATCH --mem=100000

module load cuda

python 1_Model_training_and_cell_abundance.py



