#!/usr/bin/env python
import os,sys

import kmeanconf as kc
from analysis.corelib import core
from analysis.corelib import inspect
from analysis.corelib import predict
from analysis.corelib import classifier
from analysis.corelib import optpriors
from analysis.corelib import jar
from analysis.corelib import mlsamples
from analysis.corelib import summary

from analysis.qpdlib  import tolhapdf
from analysis.qpdlib  import benchmark
from analysis.qpdlib  import pdf_pion

from analysis.obslib  import dy_pion
from analysis.obslib  import ln
from analysis.obslib  import dy_pT_pion

wdir=sys.argv[1]

##--initial processeing 
inspect.get_msr_inspected(wdir)
#predict.get_predictions(wdir)
#classifier.gen_labels(wdir,kc)
#jar.gen_jar_file(wdir,kc)
#summary.print_summary(wdir,kc)
####optpriors.gen_priors(wdir,kc,1000)
###
#####--observables
#dy_pion.plot_obs(wdir,kc)
##ln.plot_obs(wdir,kc)
##dy_pT_pion.plot_obs(wdir,kc)
###
####--pdfs
#pdf_pion.gen_xf(wdir,Q2=10.0)
#pdf_pion.plot3_xf(wdir,kc,Q2=10.0)
##pdf_pion.gen_xf(wdir,Q2=None)
##pdf_pion.plot3_xf(wdir,kc,Q2=None)
#pdf_pion.gen_moms(wdir,Q2=None)
#pdf_pion.plot1_mom1(wdir,kc,Q2=None)
#pdf_pion.print_mom1(wdir,kc,Q2=None)
#
#--pi2n
#pdf_pion.plot1_pi2n(wdir,kc)







