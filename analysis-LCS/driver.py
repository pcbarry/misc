#!/usr/bin/env python
import os,sys

import kmeanconf as kc
from analysis.corelib import core
from analysis.corelib import inspect
from analysis.corelib import predict
from analysis.corelib import classifier
#from analysis.corelib import optpriors
#from analysis.corelib import jar
#from analysis.corelib import mlsamples
from analysis.corelib import summary
#from analysis.qpdlib  import tolhapdf
#from analysis.qpdlib  import benchmark
from analysis.qpdlib  import pdf_pion

#from lib import pos_mom
#from lib import neg_mom
#from lib import dy_pion
#from lib import dy_pion2
#from lib import ln
#from lib import dy_pT_pion
#from lib import predict as pion_predict
#from lib import pdfcolor
#from lib import params
#from lib import kin
#from lib import plotchannel

try:
    basedir=sys.argv[2]
    basestep=sys.argv[3]
except:
    basedir=False
    basestep=False

wdir=sys.argv[1]

#--initial processeing 
#inspect.get_msr_inspected(wdir,FILT=[['g1 a',-1.5,'less']])
inspect.get_msr_inspected(wdir)
predict.get_predictions(wdir)
classifier.gen_labels(wdir,kc)
jar.gen_jar_file(wdir,kc)
#optpriors.gen_normal_priors(wdir,kc,1000)

#pdf_pion.gen_moms(wdir,Q2=None)
#pos_mom.sort_moms(wdir,Q2=10.0)

#--pdfs
pdf_pion.gen_xf(wdir,Q2=10.0)
#pdf_pion.plot_xf(wdir,kc,Q2=10.0)
#pdfcolor.plot_xf_expt(wdir,kc,Q2=10.0)
pdf_pion.gen_xf(wdir)
#pdf_pion.plot_xf(wdir,kc)
#pdf_pion.gen_moms(wdir,Q2=None)
#pdf_pion.plot_mom1(wdir,kc)
#pdf_pion.plot_mom_tot(wdir,kc) #--only for the NOMSR results
#pdf_pion.print_mom1(wdir,kc,Q2=None)

#--pdf ratios
#if basedir and basestep: 
#    pdf_pion.plot_xf_rat(basedir,basestep,kc,Q2=10.0)

params.plot_params(wdir)

#--pi2n
#pdf_pion.plot_pi2n(wdir,kc)

##--positive moments
#pos_mom.sort_moms(wdir,Q2=10.0)
#pos_mom.plot_xf(wdir,kc,Q2=10.0)
##pos_mom.sort_moms(wdir)
##pos_mom.plot_xf(wdir,kc)
#pos_mom.plot_mom1(wdir,kc)
#
##--negative moments
##neg_mom.sort_moms(wdir,Q2=10.0)
#
##--pos vs neg moms
##pos_mom.plot_chi2(wdir)
##neg_mom.plot_chi2(wdir)
##pion_predict.get_predictions(wdir,sign='pos')
##pion_predict.get_predictions(wdir,sign='neg')
##dy_pion2.plot_rat_pos_neg(wdir,kc)
