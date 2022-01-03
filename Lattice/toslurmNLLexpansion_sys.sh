#!/bin/csh
#SBATCH --account=f4thy
#SBATCH --nodes 1
#SBATCH --partition=production
#SBATCH --cpus-per-task 5
#SBATCH --mem=1G
#SBATCH --time=20:00:00
#SBATCH --constraint=general
#SBATCH --job-name="pITD B1B2P1P2-expansion"
#SBATCH --output=/farm_out/barryp/analysis-LCS/pITD/B1B2P1P2-expansion/out/8CW9XVS6MWVA.msr.out
#SBATCH --error=/farm_out/barryp/analysis-LCS/pITD/B1B2P1P2-expansion/out/8CW9XVS6MWVA.msr.err
/home/barryp/.conda/envs/pb3/bin/python mainNLLexpansion_sys.py 8CW9XVS6MWVA.msr B1B2P1P2-expansion pITD
