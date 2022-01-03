#!/bin/csh
#SBATCH --account=jam
#SBATCH --nodes 1
#SBATCH --partition=production
#SBATCH --cpus-per-task 5
#SBATCH --mem=1G
#SBATCH --time=20:00:00
#SBATCH --constraint=general
#SBATCH --job-name="ccLCS DMgood_stat_u_50"
#SBATCH --output=/farm_out/barryp/analysis-LCS/ccLCS/DMgood_stat_u_50/out/SUTDBX3I8AAQ.msr.npy.out
#SBATCH --error=/farm_out/barryp/analysis-LCS/ccLCS/DMgood_stat_u_50/out/SUTDBX3I8AAQ.msr.npy.err
/home/barryp/.conda/envs/pb3/bin/python  main50.py SUTDBX3I8AAQ.msr.npy DMgood_stat_u_50 ccLCS 
