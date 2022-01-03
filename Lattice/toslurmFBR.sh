#!/bin/csh
#SBATCH --account=jam
#SBATCH --nodes 1
#SBATCH --partition=production
#SBATCH --cpus-per-task 5
#SBATCH --mem=1G
#SBATCH --time=20:00:00
#SBATCH --constraint=general
#SBATCH --job-name="ccLCS DM_F_B_R"
#SBATCH --output=/farm_out/barryp/analysis-LCS/ccLCS/DM_F_B_R/out/CG592FWA502D.msr.out
#SBATCH --error=/farm_out/barryp/analysis-LCS/ccLCS/DM_F_B_R/out/CG592FWA502D.msr.err
/home/barryp/.conda/envs/pb3/bin/python  mainFBR.py CG592FWA502D.msr DM_F_B_R ccLCS 
