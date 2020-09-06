#!/bin/bash

#SBATCH -p standard

#SBATCH -J Ejecta

#SBATCH -o output_%j

#SBATCH --mem-per-cpu=8GB

#SBATCH -t 100:00:00

#SBATCH -n 1

python ejecta_search.py

python Launch.py
python Landing.py
python Panels.py
