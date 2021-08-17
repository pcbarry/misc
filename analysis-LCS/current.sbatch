#!/bin/csh
#SBATCH --account=jam
#SBATCH --nodes 1
#SBATCH --partition=production
#SBATCH --cpus-per-task 5
#SBATCH --mem=2G
#SBATCH --time=0-20:00:00
#SBATCH --constraint=general
#SBATCH --job-name="jam-6ZBKXC8RMN"
#SBATCH --output=/w/jam-sciwork18/barryp/JAM/analysis-LCS/ccLCSnew/test//out/jam-6ZBKXC8RMN.out
#SBATCH --error=/w/jam-sciwork18/barryp/JAM/analysis-LCS/ccLCSnew/test//out/jam-6ZBKXC8RMN.err

setenv FITPACK /w/jam-sciwork18/barryp/JAM/fitpack2
setenv PYTHONPATH /w/jam-sciwork18/barryp/JAM/fitpack2:/work/JAM/apps/LHAPDF/lib/python3.6/site-packages/
#cp /w/jam-sciwork18/barryp/JAM/analysis-LCS/input.py  ./
/home/barryp/.conda/envs/pb3/bin/python  /w/jam-sciwork18/barryp/JAM/fitpack2/fitlib/maxlike.py  -n 10 -v 100 input.py -msrh "mv <<fname>> /w/jam-sciwork18/barryp/JAM/analysis-LCS/ccLCSnew/test//msr/"   -p /w/jam-sciwork18/barryp/JAM/analysis-LCS/NLOpy3//msr-inspected/LCKMMJBIJ1S3.msr

