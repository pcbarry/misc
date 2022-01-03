#!/bin/csh
#SBATCH --account=jam
#SBATCH --nodes 1
#SBATCH --partition=production
#SBATCH --cpus-per-task 5
#SBATCH --mem=1G
#SBATCH --time=20:00:00
#SBATCH --constraint=general
#SBATCH --job-name="ccLCS DM_stat_u_100"
#SBATCH --output=/farm_out/barryp/analysis-LCS/ccLCS/DM_stat_u_100/out/A66AD2T361G7.msr.npy.out
#SBATCH --error=/farm_out/barryp/analysis-LCS/ccLCS/DM_stat_u_100/out/A66AD2T361G7.msr.npy.err
/home/barryp/.conda/envs/pb3/bin/python  main100.py A66AD2T361G7.msr.npy DM_stat_u_100 ccLCS 
