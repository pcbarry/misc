#!/usr/bin/env python
import sys, os
import numpy as np
import copy
from subprocess import Popen, PIPE, STDOUT
import pandas as pd
import scipy as sp
# import time
# from glob import glob

## matplotlib
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
from matplotlib.ticker import MultipleLocator, FormatStrFormatter ## for minor ticks in x label
matplotlib.rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
matplotlib.rc('text', usetex = True)
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
import matplotlib.colors as mcolors
import pylab as py

## from fitpack tools
from tools.tools     import load, save, checkdir, lprint
from tools.config    import conf, load_config

## from fitpack fitlib
from fitlib.resman import RESMAN

## from fitpack analysis
from analysis.corelib import core
from analysis.corelib import classifier

## from fitpack obslib
from obslib.jets.jet_tools import find_eta_bins

def chi2_histogram(wdir, reaction, dataset, dpi):
    ## plot histogram of relative 'chi2' for 'dataset' in 'reaction'
    load_config('%s/input.py' % wdir)
    istep = core.get_istep()
    predictions = load('%s/data/predictions-%d.dat' % (wdir, istep))
    if reaction not in predictions['reactions']:
        print('%s not in predictions' % reaction)
        return
    collaboration = predictions['reactions'][reaction][dataset]['col'][0]
    print('\nplotting relative chi2 histogram for %d in %s from %s' % (dataset, reaction, wdir))
    replicas = core.get_replicas(wdir)
    relative_chi2 = []
    for replica in replicas:
        # print replica['chi2'][reaction][dataset]
        relative_chi2.append(replica['chi2'][reaction][dataset]['chi2'] / float(replica['chi2'][reaction][dataset]['npts']))

    n_rows, n_columns = 1, 1
    figure, ax = py.subplots(n_rows, n_columns)

    ax.hist(relative_chi2, len(relative_chi2), histtype = 'bar', color = 'b', label = r'$\mathrm{%s,~%s}$' % (reaction, collaboration))
    ax.legend(fontsize = 20, loc = 'upper right')
    # ax.set_ylabel(r'\boldmath$\mathrm{n}$', size = 24)
    # ax.yaxis.set_label_coords(-0.17, 2.0)

    ax.set_xlabel(r'\boldmath$\chi^2/\mathrm{n}$', size = 24)
    # ax.xaxis.set_label_coords(0.90, -0.12)

    py.subplots_adjust(left = 0.12, bottom = 0.08, right = 0.99, top = 0.97, wspace = None, hspace = 0.2)
    py.tight_layout()
    py.savefig('%s/gallery/chi2-histogram-%s-%s-%d.png' % (wdir, reaction, collaboration, istep), dpi = dpi)

if __name__ == "__main__":
    pass
