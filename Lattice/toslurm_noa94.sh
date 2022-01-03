#!/bin/csh
#SBATCH --account=jam
#SBATCH --nodes 1
#SBATCH --partition=production
#SBATCH --cpus-per-task 5
#SBATCH --mem=1G
#SBATCH --time=20:00:00
#SBATCH --constraint=general
#SBATCH --job-name="ccLCS cut_smallz_noa94"
#SBATCH --output=/farm_out/barryp/analysis-LCS/ccLCS/cut_smallz_noa94/out/CG592FWA502D.msr.out
#SBATCH --error=/farm_out/barryp/analysis-LCS/ccLCS/cut_smallz_noa94/out/CG592FWA502D.msr.err
/home/barryp/.conda/envs/pb3/bin/python  main_noa94.py CG592FWA502D.msr cut_smallz_noa94 ccLCS 
