#!/bin/csh
#SBATCH --account=jam
#SBATCH --nodes 1
#SBATCH --partition=production
#SBATCH --cpus-per-task 5
#SBATCH --mem=2G
#SBATCH --time=10:00:00
#SBATCH --constraint=general
#SBATCH --job-name="pITD B1B2P1P2F1F2-allsys-nosmallp3"
#SBATCH --output=/w/jam-sciwork18/barryp/Notebooks-Jupyterhub/Lattice/analysis-LCS/pITD/B1B2P1P2F1F2-allsys-nosmallp3/out/SZJ4X474D7TR.msr.out
#SBATCH --error=/w/jam-sciwork18/barryp/Notebooks-Jupyterhub/Lattice/analysis-LCS/pITD/B1B2P1P2F1F2-allsys-nosmallp3/out/SZJ4X474D7TR.msr.err
/home/barryp/.conda/envs/pb3/bin/python main.py SZJ4X474D7TR.msr B1B2P1P2F1F2-allsys-nosmallp3 pITD
