#!/bin/csh
#SBATCH --account=theory
#SBATCH --nodes 1
#SBATCH --partition=theory
#SBATCH --cpus-per-task 5
#SBATCH --mem=2G
#SBATCH --time=0-20:00:00
#SBATCH --constraint=general
#SBATCH --job-name="jam-IYCYREYJ15"
#SBATCH --output=/w/general-scifs17exp/JAM/barryp/JAM/analysis-pions/resultsFinal/step13b/out/jam-IYCYREYJ15.out
#SBATCH --error=/w/general-scifs17exp/JAM/barryp/JAM/analysis-pions/resultsFinal/step13b/out/jam-IYCYREYJ15.err

setenv FITPACK /work/JAM/barryp/JAM/fitpack/
setenv PYTHONPATH /work/JAM/barryp/JAM/fitpack/:/work/JAM/apps/lhapdf2/lib/python2.7/site-packages/
#cp /w/general-scifs17exp/JAM/barryp/JAM/analysis-pions/input.py  ./
/w/general-scifs17exp/JAM//apps/anaconda2/bin/python  /work/JAM/barryp/JAM/fitpack//fitlib/maxlike.py  -n 10 -v 100 input.py -msrh "mv <<fname>> /w/general-scifs17exp/JAM/barryp/JAM/analysis-pions/resultsFinal/step13b/msr/"   -p /w/general-scifs17exp/JAM/barryp/JAM/analysis-pions/resultsFinal/step12//msr-inspected/4H39GULP6BAX.msr

