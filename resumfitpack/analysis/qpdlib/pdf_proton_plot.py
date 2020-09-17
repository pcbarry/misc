#!/usr/bin/env python
import os, sys
import subprocess
import numpy as np
import scipy as sp
import pandas as pd
import copy
import re

## matplotlib
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
matplotlib.rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
matplotlib.rc('text', usetex = True)
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
import pylab as py
# from fflib.dss.dss   import DSS
# from fflib.akk.akk   import AKK
# from fflib.hkns.hkns import HKNS
# from fflib.kkp.kkp   import KKP

## from fitpack
from tools.tools import load, save, checkdir, lprint
from tools.config import conf, load_config
from analysis.corelib import core
from analysis.corelib import classifier
# from pdfcalc  import PDFCALC

class PDF_PLOT_CORE:

    def plot_lines(self, ax, flav, c, n_replicas, lab = ''): ## plot all the PDF replicas as lines
        line_width = 10.0 / float(n_replicas)
        if line_width > 1.0: line_width = 1.0
        X = self.xf_data['X']
        n_replicas = len(self.xf_data['XF'][flav])
        for i in range(n_replicas):
            if self.cluster[i] != self.best_cluster: continue
            ax.plot(X, self.xf_data['XF'][flav][i], color = c, label = lab, linewidth = line_width)

    def plot_band_JAM(self, ax, flav, c, lab, scale = 1, distinction = {'alpha': 0.7, 'zorder': 8}):
        X = self.xf_data['X']
        n_replicas = len(self.xf_data['XF'][flav])
        Y = []
        for i in range(n_replicas):
            if self.cluster[i] != self.best_cluster: continue
            Y.append(self.xf_data['XF'][flav][i])        
        Y0 = np.mean(Y, axis = 0)
        dY = np.std(Y, axis = 0)
        alpha = distinction['alpha']
        zorder = distinction['zorder']
        return ax.fill_between(X, (Y0 - dY) * scale, (Y0 + dY) * scale, color = c, label = lab, alpha = alpha, zorder = 1)

    def plot_band(self, ax, group, name, flav, alpha = 0.2, label = '', scale = 1):
        X = group['X']
        D = group['XF'][name][flav]
        zorder = 2
        if name == 'NNPDF':
            alpha = 0.3
            zorder = 2
        y_down, y_up = scale * (D['center'] - D['difference']), scale * (D['center'] + D['difference'])
        return ax.fill_between(X, y_down, y_up, color = group['XF'][name]['color'], alpha = alpha, label = label, zorder = zorder)

    def plot_ratios(self, ax, flav, c, n_replicas, lab = ''): ## plot all the PDF replicas as separate lines
        line_width = 10.0 / float(n_replicas)
        if line_width > 1.0: line_width = 1.0
        X = self.xf_data['X']
        n_replicas = len(self.xf_data['XF'][flav])
        for i in range(n_replicas):
            if self.cluster[i] != self.best_cluster: continue
            ax.plot(X, self.xf_data['XF'][flav][i] / self.base_xf_data['mean_xf'][flav], color = c, label = lab, linewidth = line_width)

    def plot_band_steps(self, ax, step, flav, color, alpha, label = None):
        X = self.xf_data[step]['X']
        D = self.xf_data[step]['XF'][flav]
        Y = []
        n_replicas = len(self.xf_data[step]['XF'][flav])
        for i in range(n_replicas):
            if self.cluster[step][i] != self.best_cluster[step]: continue
            Y.append(self.xf_data[step]['XF'][flav][i])        
        Y0 = np.mean(Y, axis = 0)
        dY = np.std(Y, axis = 0)
        return ax.fill_between(X, (Y0 - dY), (Y0 + dY), color = color, alpha = alpha, zorder = step, label = label)

    def n_histogram(self, ax, flav, n_replicas, label): ## histogram for PDF normalization
        # for i in range(n_replicas):
        #     if self.cluster[i] != self.best_cluster: continue
        if flav == 'uv':
            ax.hist(self.norm_data['first'][flav + '1'], n_replicas, histtype = 'bar', color = 'r', label = label.rsplit('$', 1)[0] + '^{(1)}$')
            ax.hist(self.norm_data['second'][flav + '2'], n_replicas, histtype = 'bar', color = 'b', label = label.rsplit('$', 1)[0] + '^{(2)}$')
        else:
            ax.hist(self.norm_data['first'][flav + '1'], n_replicas, histtype = 'bar', color = 'r', label = label)
            ax.hist(self.norm_data['second'][flav + '2'], n_replicas, histtype = 'bar', color = 'b')

class PLOTS:

    def __init__(self, task, wdir,wdir2,kc, Q2 = None, dpi = 200):

        if  task == 1:
            self.msg = 'mplots.line_figure'
            self.line_figure(wdir,wdir2,kc, Q2, dpi)
            # self.loop_over_steps(wdir,'simple',last)

        if  task == 2:
            self.msg = 'mplots.groups_comparison_figure'
            self.groups_comparison_figure(wdir, Q2, dpi)
            # self.loop_over_steps(wdir,'simple',last)

        if  task == 3:
            self.msg = 'mplots.ratio_figure'
            self.ratio_figure(wdir, Q2, dpi)
            # self.loop_over_steps(wdir,'simple',last)
            ## 'wdir' for 'self.ratio_figure' is a dictionary that contains keys 'current' and 'base'
            ## 'self.ratio_figure' will plot ratio of 'current' 'xf' values over 'base' 'xf' values
            ## wdir = {'base': './results1/step04/', 'current': './results1/step06/'}
            ## input file of current step will be loaded

        if  task == 4:
            self.msg = 'mplots.steps_comparison_figure'
            self.steps_comparison_figure(wdir, Q2, dpi)
            # self.loop_over_steps(wdir,'simple',last)
            ## 'wdir' for 'self.ratio_figure' is a dictionary that contains integer keys
            ## the keys have to be the number of the steps
            ## 'xf' of larger key values will be plotted on top of that of smaller key values
            ## 'xf' of larger key values will also be plotted with less opacity
            ## wdir = {2: './results1/step02/', 3: './results1/step03/', 4: './results1/step04/'}
            ## input file of largest key value will be loaded

        if  task == 10:
            self.msg = 'mplots.normalization_histogram'
            self.histogram_figure(wdir, Q2, dpi)
            # self.loop_over_steps(wdir,'simple',last)

    def plot_lines(self, ax, flav, c, n_replicas, lab = ''): ## plot all the PDF replicas as separate lines
        line_width = 10.0 / float(n_replicas)
        if line_width > 1.0: line_width = 1.0
        X = self.xf_data['X']
        #n_replicas = len(self.xf_data['XF'][flav])
        for i in range(n_replicas):
            #if self.cluster[i] != self.best_cluster: continue
            ax.plot(X, self.xf_data['XF'][flav][i], color = self.colors[self.cluster[i]], label = lab,alpha=0.2)#, linewidth = line_width)

    def plot_band_JAM(self, ax, flav, c, lab, scale = 1):
        X = self.xf_data['X']
        n_replicas = len(self.xf_data['XF'][flav])
        Y = []
        for i in range(n_replicas):
            if self.cluster[i] != self.best_cluster: continue
            Y.append(self.xf_data['XF'][flav][i])        
        Y0 = np.mean(Y, axis = 0)
        dY = np.std(Y, axis = 0)
        return ax.fill_between(X, (Y0 - dY) * scale, (Y0 + dY) * scale, color = c, label = lab, alpha = 0.7, zorder = 8)

    def plot_band(self, ax, group, name, flav, alpha = 0.2, label = '', scale = 1):
        X = group['X']
        D = group['XF'][name][flav]
        zorder = 2
        if name == 'NNPDF':
            alpha = 0.3
            zorder = 2
        return ax.fill_between(X, scale * (D['center'] - D['difference']), scale * (D['center'] + D['difference']), color = group['XF'][name]['color'], alpha = alpha, label = label, zorder = zorder)

    def plot_ratios(self, ax, flav, c, n_replicas, lab = ''): ## plot all the PDF replicas as separate lines
        line_width = 10.0 / float(n_replicas)
        if line_width > 1.0: line_width = 1.0
        X = self.xf_data['X']
        n_replicas = len(self.xf_data['XF'][flav])
        for i in range(n_replicas):
            if self.cluster[i] != self.best_cluster: continue
            ax.plot(X, self.xf_data['XF'][flav][i] / self.base_xf_data['mean_xf'][flav], color = c, label = lab, linewidth = line_width)

    def plot_band_steps(self, ax, step, flav, color, alpha, label = None):
        X = self.xf_data[step]['X']
        D = self.xf_data[step]['XF'][flav]
        Y = []
        n_replicas = len(self.xf_data[step]['XF'][flav])
        for i in range(n_replicas):
            if self.cluster[step][i] != self.best_cluster[step]: continue
            Y.append(self.xf_data[step]['XF'][flav][i])        
        Y0 = np.mean(Y, axis = 0)
        dY = np.std(Y, axis = 0)
        return ax.fill_between(X, (Y0 - dY), (Y0 + dY), color = color, alpha = alpha, zorder = step, label = label)

    def n_histogram(self, ax, flav, n_replicas, label): ## plot all the PDF normalization parameters
        # for i in range(n_replicas):
        #     if self.cluster[i] != self.best_cluster: continue
        if flav == 'uv':
            ax.hist(self.norm_data['first'][flav + '1'], n_replicas, histtype = 'bar', color = 'r', label = label.rsplit('$', 1)[0] + '^{(1)}$')
            ax.hist(self.norm_data['second'][flav + '2'], n_replicas, histtype = 'bar', color = 'b', label = label.rsplit('$', 1)[0] + '^{(2)}$')
        else:
            ax.hist(self.norm_data['first'][flav + '1'], n_replicas, histtype = 'bar', color = 'r', label = label)
            ax.hist(self.norm_data['second'][flav + '2'], n_replicas, histtype = 'bar', color = 'b')

    def line_figure(self, wdir,wdir2,kc, Q2, dpi):
        load_config('%s/input.py' % wdir)
        istep = core.get_istep()
        #labels = load('%s/data/labels-%d.dat' % (wdir, istep))
        #self.cluster = labels['cluster']
        #n_replicas = len(self.cluster)
        #self.best_cluster = 1

        nrows, ncols = 4, 2
        fig = py.figure(figsize = (ncols * 5.0, nrows * 3.0))
        AX = [py.subplot(nrows, ncols, cnt+1) for cnt in range(8) if cnt != 7]

        ## PDF
        # self.xf_data = load('%s/data/pdf-%d.dat' % (wdir, istep))
        if Q2 == None: Q2 = conf['Q20']
        if Q2 == conf['Q20']:
            self.xf_data = load('%s/data/pdf-%d.dat' % (wdir, istep))
        else:
            self.xf_data = load('%s/data/pdf-%d-%f.dat' % (wdir, istep, Q2))

        load_config('%s/input.py'%wdir2)
        istep2=core.get_istep()
        self.cluster,self.colors,self.nc,self.cluster_order = classifier.get_clusters(wdir2,istep2,kc)
        self.best_cluster=self.cluster_order[0]
        n_replicas=len(self.cluster)

        print '\nplotting PDF line figure from %s' % wdir
        self.plot_lines(AX[0], 'uv', 'r', n_replicas)
        self.plot_lines(AX[0], 'dv', 'b', n_replicas)
        self.plot_lines(AX[1], 'd/u' , 'r', n_replicas)
        self.plot_lines(AX[2], 'db+ub', 'r', n_replicas)
        self.plot_lines(AX[3], 'db-ub', 'r', n_replicas)
        self.plot_lines(AX[4], 's+sb', 'r', n_replicas)
        self.plot_lines(AX[5], 'rs', 'r', n_replicas)
        self.plot_lines(AX[6], 'g', 'r', n_replicas)

        AX[0].set_ylim(0.0, 0.85)
        AX[0].set_yticks([0.0, 0.2, 0.4, 0.6, 0.8])
        AX[0].set_yticklabels([r'$0$', r'$0.2$', r'$0.4$', r'$0.6$', r'$0.8$'])
        AX[0].text(0.4, 0.81, r'\boldmath$xu_v$', color = 'k', transform = AX[0].transAxes, size = 28)
        AX[0].text(0.4, 0.07, r'\boldmath$xd_v$', color = 'k', transform = AX[0].transAxes, size = 28)
        AX[1].set_ylim(0.0, 1.0)
        AX[1].set_yticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
        AX[1].set_yticklabels([r'$0$', r'', r'$0.2$', r'', r'$0.4$', r'', r'$0.6$', r'', r'$0.8$'])
        AX[1].text(0.07, 0.13, r'\boldmath$d/u$', color = 'k', transform = AX[1].transAxes, size = 28)
        AX[2].set_ylim(0.0, 0.5)
        AX[2].set_yticks([0.0, 0.1, 0.2, 0.3, 0.4])
        AX[2].set_yticklabels([r'$0$', r'', r'$0.2$', r'', r'$0.4$'])
        AX[2].text(0.07, 0.1, r'\boldmath$x(\bar{d}\!+\!\bar{u})$', color = 'k', transform = AX[2].transAxes, size = 28)
        AX[3].set_ylim(-0.05, 0.13)
        AX[3].set_yticks([-0.05, 0.0, 0.05, 0.1])
        AX[3].set_yticklabels([r'-$0.05$', r'$0$', r'$0.05$', r'$0.1$'])
        AX[3].text(0.07, 0.81, r'\boldmath$x(\bar{d}\!-\!\bar{u})$', color = 'k', transform = AX[3].transAxes, size = 28)
        AX[4].set_ylim(0.0, 0.48)
        AX[4].set_yticks([0.0, 0.1, 0.2, 0.3, 0.4])
        AX[4].set_yticklabels([r'$0$', r'', r'$0.2$', r'', r'$0.4$'])
        AX[4].text(0.6, 0.8, r'\boldmath$x(s\!+\!\bar{s})$', color = 'k', transform = AX[4].transAxes, size = 28)
        AX[5].set_ylim(0.0, 1.6)
        AX[5].text(0.05, 0.78, r'\boldmath$R_s$', color = 'k', transform = AX[5].transAxes, size = 30)
        AX[5].set_yticks([0.0, 0.5, 1.0, 1.5])
        AX[5].set_yticklabels([r'$0$', r'$0.5$', r'$1$', r'$1.5$'])
        AX[6].set_ylim(0.0, 2.5)
        AX[6].text(0.8, 0.8, r'\boldmath$xg$', color = 'k', transform = AX[6].transAxes, size = 31)
        AX[6].set_yticks([0.0, 0.5, 1.0, 1.5, 2.0])

        for ax in AX:
            ax.semilogx()
            ax.set_xlim(8e-3, 9e-1)
            ax.set_xticks([1e-2, 1e-1])
            ax.xaxis.set_label_coords(0.95, -0.05)
            ax.tick_params(axis = 'both', which = 'both', right = True, top = True, direction = 'in', labelsize = 20)
            ax.set_xticklabels([r'', r''])
            ax.legend(frameon = 0, loc = 'best')

        AX[0].set_ylabel(r'$xf(x)$', size = 20)
        AX[5].set_xlabel(r'\boldmath$x$', size = 30)
        AX[5].set_xticklabels([r'$0.01$', r'$0.1$'])
        AX[6].set_xlabel(r'\boldmath$x$', size = 30)
        AX[6].set_xticklabels([r'$0.01$', r'$0.1$'])

        py.subplots_adjust(left = 0.12, bottom = 0.08, right = 0.99, top = 0.97, wspace = None, hspace = 0.2)
        py.tight_layout()
        checkdir('%s/gallery'%wdir)
        if Q2 == conf['Q20']:
            py.savefig('%s/gallery/pdf-lines-%d.png' % (wdir, istep), dpi = dpi)
        else:
            py.savefig('%s/gallery/pdf-lines-%d-%f.png' % (wdir, istep, Q2), dpi = dpi)

    def groups_comparison_figure(self, wdir, Q2, dpi):
        load_config('%s/input.py' % wdir)
        istep = core.get_istep()
        labels = load('%s/data/labels-%d.dat' % (wdir, istep))
        self.cluster = labels['cluster']
        n_replicas = len(self.cluster)
        self.best_cluster = 0

        if Q2 == None: Q2 = conf['Q20']
        groups = load('%s/analysis/qpdlib/lhapdf-%f.dat' % (os.environ['FITPACK'], Q2))

        nrows, ncols = 4, 2
        fig = py.figure(figsize = (ncols * 5.0, nrows * 3.0))
        AX = [py.subplot(nrows, ncols, cnt+1) for cnt in range(8) if cnt != 7]

        ## PDF
        if Q2 == conf['Q20']:
            self.xf_data = load('%s/data/pdf-%d.dat' % (wdir, istep))
        else:
            self.xf_data = load('%s/data/pdf-%d-%f.dat' % (wdir, istep, Q2))

        print '\nplotting PDF band figure with following groups at Q2 = %f' % Q2
        print groups['groups']

        ax = AX[0]
        self.plot_band_JAM(ax, 'uv', 'r', r'$\rm JAM$')
        for _ in groups['XF']:
            self.plot_band(ax, groups, _, 'uv', alpha = 0.2, label = r'$\rm %s$' % _)
        self.plot_band_JAM(ax, 'dv', 'r', r'$\rm JAM$')
        for _ in groups['XF']:
            self.plot_band(ax, groups, _, 'dv', alpha = 0.2, label = r'$\rm %s$' % _)
        ax.text(0.4, 0.81, r'\boldmath$xu_v$', color = 'k', transform = ax.transAxes, size = 28)
        ax.text(0.6, 0.1, r'\boldmath$xd_v$', color = 'k', transform = ax.transAxes, size = 28)
        ax.set_ylim(0.0, 0.85)
        ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8])
        ax.set_yticklabels([r'$0$', r'$0.2$', r'$0.4$', r'$0.6$', r'$0.8$'])

        ax = AX[1]
        self.plot_band_JAM(ax, 'd/u', 'r', r'$\rm JAM$')
        for _ in groups['XF']:
            self.plot_band(ax, groups, _, 'd/u', alpha = 0.2, label = r'$\rm %s$' % _)
        ax.text(0.07, 0.13, r'\boldmath$d/u$', color = 'k', transform = ax.transAxes, size = 30)
        ax.set_ylim(0.0, 0.9)
        ax.set_yticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
        ax.set_yticklabels([r'$0$', r'', r'$0.2$', r'', r'$0.4$', r'', r'$0.6$', r'', r'$0.8$'])

        ax = AX[2]
        self.plot_band_JAM(ax, 'db+ub','r', r'$\rm JAM$')
        for _ in groups['XF']:
            self.plot_band(ax, groups, _, 'db+ub', alpha = 0.2, label = r'$\rm %s$' % _)
        ax.text(0.07, 0.13, r'\boldmath$x(\bar{d}\!+\!\bar{u})$', color = 'k', transform = ax.transAxes, size = 28)
        ax.set_ylim(0.0, 0.5)
        ax.set_yticks([0.0, 0.1, 0.2, 0.3, 0.4])
        ax.set_yticklabels([r'$0$', r'', r'$0.2$', r'', r'$0.4$'])

        ax = AX[3]
        self.plot_band_JAM(ax, 'db-ub', 'r', r'$\rm JAM$')
        for _ in groups['XF']:
            self.plot_band(ax, groups, _, 'db-ub', alpha = 0.2, label = r'$\rm %s$' % _)
        ax.text(0.07, 0.8, r'\boldmath$x(\bar{d}\!-\!\bar{u})$', color = 'k', transform = ax.transAxes,size=28)
        ax.axhline(y = 0.0, c = 'k', ls = '-', lw = 1, alpha = 0.5)
        ax.set_ylim(-0.07, 0.13)
        ax.set_yticks([-0.05, 0.0, 0.05, 0.1])
        ax.set_yticklabels([r'-$0.05$', r'$0$', r'$0.05$', r'$0.1$'])

        ax = AX[4]
        self.plot_band_JAM(ax, 's+sb', 'r', r'$\rm JAM$')
        for _ in groups['XF']:
            self.plot_band(ax, groups, _, 's+sb', alpha = 0.2, label = r'$\rm %s$' % _)
        ax.text(0.6, 0.8, r'\boldmath$x(s\!+\!\bar{s})$', color = 'k', transform = ax.transAxes, size = 28)
        ax.set_ylim(0.0, 0.48)
        ax.set_yticks([0.0, 0.1, 0.2, 0.3, 0.4])
        ax.set_yticklabels([r'$0$', r'', r'$0.2$', r'', r'$0.4$'])

        ax = AX[5]
        self.plot_band_JAM(ax, 'rs', 'r', r'$\rm JAM$')
        for _ in groups['XF']:
            self.plot_band(ax, groups, _, 'rs', alpha = 0.2, label = r'$\rm %s$' % _)
        ax.text(0.05, 0.78, r'\boldmath$R_s$', color = 'k', transform = ax.transAxes, size = 30)
        ax.set_ylim(0.0, 1.6)
        ax.set_yticks([0.0, 0.5, 1.0, 1.5])
        ax.set_yticklabels([r'$0$', r'$0.5$', r'$1$', r'$1.5$'])

        ax = AX[6]
        self.plot_band_JAM(ax, 'g', 'r', r'$\rm JAM~Jet$')
        for _ in groups['XF']:
            if _ == 'CJ15nlo': label=r'$\rm CJ15$'
            elif _ == 'NNPDF31_nlo_as_0118': label = r'$\rm NNPDF3.1$'
            elif _ == 'MMHT2014nlo68cl': label = r'$\rm MMHT14$'
            elif _ == 'ABMP16_3_nlo': label = r'$\rm ABMP16$'
            elif _ == 'CSKK_nnlo_EIG': label = r'$\rm CSSK$'
            elif _ == 'JAM19PDF_proton_nlo': label = r'$\rm JAM19$'
            else: label = r'$\rm %s$' % _
            self.plot_band(ax, groups, _, 'g', alpha = 0.2, label = label)
        ax.text(0.8, 0.8, r'\boldmath$xg$', color = 'k', transform = ax.transAxes, size = 31)
        ax.set_ylim(0.0, 3.7)
        # ax.set_yticks([0.0, 0.5, 1.0, 1.5, 2.0])
        ax.set_yticks([0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5])
        ax.legend(frameon = 0, loc = 'best', bbox_to_anchor = [1.87, 0.9], fontsize = 15)

        for ax in AX:
            ax.semilogx()
            ax.tick_params(axis = 'both', which = 'both', right = True, top = True, direction = 'in', labelsize = 20)
            ax.xaxis.set_label_coords(0.7, -0.05)
            ax.set_xlim(8e-3, 9e-1)
            ax.set_xticks([1e-2, 1e-1])
            if ax == AX[5] or ax == AX[6]:
                ax.set_xlabel(r'\boldmath$x$', size = 30)   
                ax.set_xticklabels([r'$0.01$', r'$0.1$'])
                ax.text(0.84, -0.095, r'$0.5$', color = 'k', transform = ax.transAxes, size = 20)
            else:
                ax.set_xticklabels([r'', r''])

        py.tight_layout()
        py.subplots_adjust(left=0.05, bottom=0.04, right=0.99, top=0.99, wspace=0.17, hspace=0.04)
        if Q2 == conf['Q20']:
            py.savefig('%s/gallery/pdf-groups-%d.png' % (wdir, istep), dpi = dpi)
        else:
            py.savefig('%s/gallery/pdf-groups-%d-%f.png' % (wdir, istep, Q2), dpi = dpi)

    def ratio_figure(self, wdirs, Q2, dpi):
        load_config('%s/input.py' % wdirs['current'])
        istep = core.get_istep()
        labels = load('%s/data/labels-%d.dat' % (wdirs['current'], istep))
        self.cluster = labels['cluster']
        n_replicas = len(self.cluster)
        self.best_cluster = 0
        base_step = int(re.findall(r'\d+', wdirs['base'].split('/')[2])[0])

        nrows, ncols = 4, 2
        fig = py.figure(figsize = (ncols * 5.0, nrows * 3.0))
        AX = [py.subplot(nrows, ncols, cnt+1) for cnt in range(8) if cnt != 7]

        ## PDF
        if Q2 == None: Q2 = conf['Q20']
        if Q2 == conf['Q20']:
            self.xf_data = load('%s/data/pdf-%d.dat' % (wdirs['current'], istep))
            self.base_xf_data = load('%s/data/pdf-%d.dat' % (wdirs['base'], base_step))
        else:
            self.xf_data = load('%s/data/pdf-%d-%f.dat' % (wdirs['current'], istep, Q2))
            self.base_xf_data = load('%s/data/pdf-%d-%f.dat' % (wdirs['base'], base_step, Q2))

        self.base_xf_data['mean_xf'] = {}
        for flavor in self.base_xf_data['XF']:
            self.base_xf_data['mean_xf'][flavor] = np.mean(self.base_xf_data['XF'][flavor], axis = 0)

        print '\nplotting PDF ratio figure of %s over %s at Q2 = %f' % (wdirs['current'], wdirs['base'], Q2)
        self.plot_ratios(AX[0], 'uv', 'r', n_replicas)
        self.plot_ratios(AX[0], 'dv', 'b', n_replicas)
        self.plot_ratios(AX[1], 'd/u' , 'r', n_replicas)
        self.plot_ratios(AX[2], 'db+ub', 'r', n_replicas)
        self.plot_ratios(AX[3], 'db-ub', 'r', n_replicas)
        self.plot_ratios(AX[4], 's+sb', 'r', n_replicas)
        self.plot_ratios(AX[5], 'rs', 'r', n_replicas)
        self.plot_ratios(AX[6], 'g', 'r', n_replicas)

        AX[0].text(0.4, 0.81, r'\boldmath$xu_v$', color = 'r', transform = AX[0].transAxes, size = 28)
        AX[0].text(0.4, 0.07, r'\boldmath$xd_v$', color = 'b', transform = AX[0].transAxes, size = 28)
        AX[1].text(0.07, 0.13, r'\boldmath$d/u$', color = 'k', transform = AX[1].transAxes, size = 28)
        AX[2].text(0.07, 0.1, r'\boldmath$x(\bar{d}\!+\!\bar{u})$', color = 'k', transform = AX[2].transAxes, size = 28)
        AX[3].text(0.07, 0.81, r'\boldmath$x(\bar{d}\!-\!\bar{u})$', color = 'k', transform = AX[3].transAxes, size = 28)
        AX[4].text(0.07, 0.8, r'\boldmath$x(s\!+\!\bar{s})$', color = 'k', transform = AX[4].transAxes, size = 28)
        AX[5].text(0.05, 0.78, r'\boldmath$R_s$', color = 'k', transform = AX[5].transAxes, size = 30)
        AX[6].text(0.07, 0.8, r'\boldmath$xg$', color = 'k', transform = AX[6].transAxes, size = 31)

        for ax in AX:
            ax.semilogx()
            ax.set_xlim(8e-3, 9e-1)
            ax.set_ylim(0.8, 1.2)
            ax.set_yticks([0.9, 1.0, 1.1])
            ax.set_yticklabels([r'$0.9$', r'$1.0$', r'$1.1$'])
            ax.set_xticks([1e-2, 1e-1])
            ax.xaxis.set_label_coords(0.95, -0.05)
            ax.tick_params(axis = 'both', which = 'both', right = True, top = True, direction = 'in', labelsize = 20)
            ax.set_xticklabels([r'', r''])
            ax.axhline(1.0, color = 'k', linestyle = ':')
            ax.legend(frameon = 0, loc = 'best')

        AX[0].set_ylabel(r'$xf(x)$', size = 20)
        AX[5].set_xlabel(r'\boldmath$x$', size = 30)
        AX[5].set_xticklabels([r'$0.01$', r'$0.1$'])
        AX[6].set_xlabel(r'\boldmath$x$', size = 30)
        AX[6].set_xticklabels([r'$0.01$', r'$0.1$'])

        py.subplots_adjust(left = 0.12, bottom = 0.08, right = 0.99, top = 0.97, wspace = None, hspace = 0.2)
        py.tight_layout()
        if Q2 == conf['Q20']:
            py.savefig('%s/gallery/pdf-ratio-%d.png' % (wdirs['current'], istep), dpi = dpi)
        else:
            py.savefig('%s/gallery/pdf-ratio-%d-%f.png' % (wdirs['current'], istep, Q2), dpi = dpi)

    def steps_comparison_figure(self, wdirs, Q2, dpi):
        step_keys = sorted(wdirs.keys())
        load_config('%s/input.py' % wdirs[step_keys[-1]])
        istep = core.get_istep()

        nrows, ncols = 4, 2
        fig = py.figure(figsize = (ncols * 5.0, nrows * 3.0))
        AX = [py.subplot(nrows, ncols, cnt+1) for cnt in range(8) if cnt != 7]
        colors = ['coral', 'grey', 'wheat', 'blueviolet', 'brown', 'olivedrab', 'deepskyblue']
        # colors = ['k', 'r', 'g', 'Yellow', 'b']

        ## PDF
        self.xf_data = {}
        self.cluster = {}
        self.best_cluster = {}
        labels = {}
        if Q2 == None: Q2 = conf['Q20']
        if Q2 == conf['Q20']:
            for step_key in step_keys:
                self.xf_data[step_key] = load('%s/data/pdf-%d.dat' % (wdirs[step_key], step_key))
                labels[step_key] = load('%s/data/labels-%d.dat' % (wdirs[step_key], step_key))
                self.cluster[step_key] = labels[step_key]['cluster']
                self.best_cluster[step_key] = 0
        else:
            for step_key in step_keys:
                self.xf_data[step_key] = load('%s/data/pdf-%d-%f.dat' % (wdirs[step_key], step_key, Q2))
                labels[step_key] = load('%s/data/labels-%d.dat' % (wdirs[step_key], step_key))
                self.cluster[step_key] = labels[step_key]['cluster']
                self.best_cluster[step_key] = 0

        print '\nplotting PDF band figure with following steps at Q2 = %f' % Q2
        print [wdirs[_] for _ in wdirs]

        ax = AX[0]
        for i in range(len(step_keys)):
            step = step_keys[i]
            # alpha = (1.0 / len(step_keys)) * (len(step_keys) - i)
            alpha = 0.5
            color = colors[i]
            self.plot_band_steps(ax, step, 'uv', color, alpha)
            self.plot_band_steps(ax, step, 'dv', color, alpha)
        ax.text(0.4, 0.81, r'\boldmath$xu_v$', color = 'k', transform = ax.transAxes, size = 28)
        ax.text(0.6, 0.1, r'\boldmath$xd_v$', color = 'k', transform = ax.transAxes, size = 28)
        ax.set_ylim(0.0, 0.85)
        ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8])
        ax.set_yticklabels([r'$0$', r'$0.2$', r'$0.4$', r'$0.6$', r'$0.8$'])

        ax = AX[1]
        for i in range(len(step_keys)):
            step = step_keys[i]
            # alpha = (1.0 / len(step_keys)) * (len(step_keys) - i)
            alpha = 0.5
            color = colors[i]
            self.plot_band_steps(ax, step, 'd/u', color, alpha)
        ax.text(0.07, 0.13, r'\boldmath$d/u$', color = 'k', transform = ax.transAxes, size = 30)
        ax.set_ylim(0.0, 0.9)
        ax.set_yticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
        ax.set_yticklabels([r'$0$', r'', r'$0.2$', r'', r'$0.4$', r'', r'$0.6$', r'', r'$0.8$'])

        ax = AX[2]
        for i in range(len(step_keys)):
            step = step_keys[i]
            # alpha = (1.0 / len(step_keys)) * (len(step_keys) - i)
            alpha = 0.5
            color = colors[i]
            self.plot_band_steps(ax, step, 'db+ub', color, alpha)
        ax.text(0.07, 0.13, r'\boldmath$x(\bar{d}\!+\!\bar{u})$', color = 'k', transform = ax.transAxes, size = 28)
        ax.set_ylim(0.0, 0.5)
        ax.set_yticks([0.0, 0.1, 0.2, 0.3, 0.4])
        ax.set_yticklabels([r'$0$', r'', r'$0.2$', r'', r'$0.4$'])

        ax = AX[3]
        for i in range(len(step_keys)):
            step = step_keys[i]
            # alpha = (1.0 / len(step_keys)) * (len(step_keys) - i)
            alpha = 0.5
            color = colors[i]
            self.plot_band_steps(ax, step, 'db-ub', color, alpha)
        ax.text(0.07, 0.8, r'\boldmath$x(\bar{d}\!-\!\bar{u})$', color = 'k', transform = ax.transAxes,size=28)
        ax.axhline(y = 0.0, c = 'k', ls = '-', lw = 1, alpha = 0.5)
        ax.set_ylim(-0.07, 0.13)
        ax.set_yticks([-0.05, 0.0, 0.05, 0.1])
        ax.set_yticklabels([r'-$0.05$', r'$0$', r'$0.05$', r'$0.1$'])

        ax = AX[4]
        for i in range(len(step_keys)):
            step = step_keys[i]
            # alpha = (1.0 / len(step_keys)) * (len(step_keys) - i)
            alpha = 0.5
            color = colors[i]
            self.plot_band_steps(ax, step, 's+sb', color, alpha)
        ax.text(0.6, 0.8, r'\boldmath$x(s\!+\!\bar{s})$', color = 'k', transform = ax.transAxes, size = 28)
        ax.set_ylim(0.0, 0.48)
        ax.set_yticks([0.0, 0.1, 0.2, 0.3, 0.4])
        ax.set_yticklabels([r'$0$', r'', r'$0.2$', r'', r'$0.4$'])

        ax = AX[5]
        for i in range(len(step_keys)):
            step = step_keys[i]
            # alpha = (1.0 / len(step_keys)) * (len(step_keys) - i)
            alpha = 0.5
            color = colors[i]
            self.plot_band_steps(ax, step, 'rs', color, alpha)
        ax.text(0.05, 0.78, r'\boldmath$R_s$', color = 'k', transform = ax.transAxes, size = 30)
        ax.set_ylim(0.0, 1.6)
        ax.set_yticks([0.0, 0.5, 1.0, 1.5])
        ax.set_yticklabels([r'$0$', r'$0.5$', r'$1$', r'$1.5$'])

        ax = AX[6]
        for i in range(len(step_keys)):
            step = step_keys[i]
            # alpha = (1.0 / len(step_keys)) * (len(step_keys) - i)
            alpha = 0.5
            color = colors[i]
            label = r'$\mathrm{step~%02d}$' % step
            self.plot_band_steps(ax, step, 'g', color, alpha, label)
        ax.text(0.8, 0.8, r'\boldmath$xg$', color = 'k', transform = ax.transAxes, size = 31)
        ax.set_ylim(0.0, 2.5)
        ax.set_yticks([0.0, 0.5, 1.0, 1.5, 2.0])
        ax.legend(frameon = 0, loc = 'best', bbox_to_anchor = [1.87, 0.9], fontsize = 20)

        for ax in AX:
            ax.semilogx()
            ax.tick_params(axis = 'both', which = 'both', right = True, top = True, direction = 'in', labelsize = 20)
            ax.xaxis.set_label_coords(0.7, -0.05)
            ax.set_xlim(8e-3, 9e-1)
            ax.set_xticks([1e-2, 1e-1])
            if ax == AX[5] or ax == AX[6]:
                ax.set_xlabel(r'\boldmath$x$', size = 30)   
                ax.set_xticklabels([r'$0.01$', r'$0.1$'])
                ax.text(0.84, -0.095, r'$0.5$', color = 'k', transform = ax.transAxes, size = 20)
            else:
                ax.set_xticklabels([r'', r''])

        py.tight_layout()
        py.subplots_adjust(left=0.05, bottom=0.04, right=0.99, top=0.99, wspace=0.17, hspace=0.04)
        if Q2 == conf['Q20']:
            py.savefig('%s/gallery/pdf-steps-%d.png' % (wdirs[step_keys[-1]], istep), dpi = dpi)
        else:
            py.savefig('%s/gallery/pdf-steps-%d-%f.png' % (wdirs[step_keys[-1]], istep, Q2), dpi = dpi)

    def histogram_figure(self, wdir, Q2, dpi):
        ## plot histogram for normalizaton of all flavors except sea
        ## because first and second shape sea are fixed to be the same
        load_config('%s/input.py' % wdir)
        istep = core.get_istep()
        # jar = load('%s/data/jar-%d.dat' % (wdir, istep))
        labels = load('%s/data/labels-%d.dat' % (wdir, istep))
        self.cluster = labels['cluster']
        n_replicas = len(self.cluster)
        self.best_cluster = 0

        if Q2 == None: Q2 = conf['Q20']
        if Q2 == conf['Q20']:
            self.norm_data = load('%s/data/pdf-normalization-%d.dat' % (wdir, istep))
        else:
            self.norm_data= load('%s/data/pdf-normalization-%d-%f.dat' % (wdir, istep, Q2))

        # replicas = core.get_replicas(wdir)
        # core.mod_conf(istep, replicas[0]) ## set conf as specified in istep

        # if 'pdf' not in conf['steps'][istep]['active distributions']:
        #     print('pdf-proton not in active distribution')
        #     return

        # resman = RESMAN(nworkers = 1, parallel = False, datasets = False)
        # parman = resman.parman
        # parman.order = replicas[0]['order'][istep] ## make sure 'parman' uses the same order as all the replicas do
        # # print parman.order

        # pdf = conf['pdf']

        # core.mod_conf(istep, replicas[i_replica - 1])
        # parman.set_new_params(replicas[i_replica - 1]['params'][istep], initial = True)
        # for name, value in pdf.params.iteritems():
        #     print '%7s: %.5e, %.5e, %.5e, %.5e, %.5e' % (name, value[0], value[1], value[2], value[3], value[4])

        nrows, ncols = 4, 2
        fig = py.figure(figsize = (ncols * 5.0, nrows * 3.0))
        # AX = [py.subplot(nrows, ncols, cnt) for cnt in range(1, 9)]
        AX = [py.subplot(nrows, ncols, cnt+1) for cnt in range(8) if cnt != 7]

        ## PDF
        print '\nplotting PDF normalization histogram from %s' % wdir
        self.n_histogram(AX[0], 'uv', n_replicas, r'\boldmath$u_v$')
        self.n_histogram(AX[1], 'dv', n_replicas, r'\boldmath$d_v$')
        self.n_histogram(AX[2], 'ub' , n_replicas, r'\boldmath$u_b$')
        self.n_histogram(AX[3], 'db', n_replicas, r'\boldmath$d_b$')
        self.n_histogram(AX[4], 's', n_replicas, r'\boldmath$s$')
        self.n_histogram(AX[5], 'sb', n_replicas, r'\boldmath$s_b$')
        self.n_histogram(AX[6], 'g', n_replicas, r'\boldmath$g$')
        # self.n_histogram(AX[7], 'sea', n_replicas, r'\boldmath$sea$')

        for ax in AX:
            # ax.semilogx()
            # ax.set_xlim(8e-3, 9e-1)
            # ax.set_xticks([1e-2, 1e-1])
            # ax.set_xlabel(r'\boldmath$x$',size=15)
            # ax.xaxis.set_label_coords(0.95, -0.05)
            # ax.tick_params(axis = 'both', which = 'both', right = True, top = True, direction = 'in', labelsize = 20)
            # ax.set_xticklabels([r'', r''])
            ax.legend(frameon = 0, loc = 'best', fontsize = 17)

        ## log scale x axis has to be disabled for the time being
        ## because we need to try negative values for normalization of second shape
        # AX[2].semilogx()
        # AX[3].semilogx()
        # AX[4].semilogx()
        # AX[5].semilogx()

        py.subplots_adjust(left = 0.12, bottom = 0.08, right = 0.99, top = 0.97, wspace = None, hspace = 0.2)
        py.tight_layout()
        if Q2 == conf['Q20']:
            py.savefig('%s/gallery/pdf-normalization-histogram-%d.png' % (wdir, istep), dpi = dpi)
        else:
            py.savefig('%s/gallery/pdf-normalization-histogram-%d-%f.png' % (wdir, istep, Q2), dpi = dpi)

class NO_STRANGE(PDF_PLOT_CORE):

    def __init__(self, task, wdir, Q2 = None, dpi = 200):

        if  task == 1:
            self.msg = 'mplots.line_figure'
            self.line_figure(wdir, Q2, dpi)
            # self.loop_over_steps(wdir,'simple',last)

        if  task == 2:
            self.msg = 'mplots.groups_comparison_figure'
            self.groups_comparison_figure(wdir, Q2, dpi)
            # self.loop_over_steps(wdir,'simple',last)

        if  task == 3:
            self.msg = 'mplots.ratio_figure'
            self.ratio_figure(wdir, Q2, dpi)
            # self.loop_over_steps(wdir,'simple',last)
            ## 'wdir' for 'self.ratio_figure' is a dictionary that contains keys 'current' and 'base'
            ## 'self.ratio_figure' will plot ratio of 'current' 'xf' values over 'base' 'xf' values
            ## wdir = {'base': './results1/step04/', 'current': './results1/step06/'}
            ## input file of current step will be loaded

        if  task == 4:
            self.msg = 'mplots.steps_comparison_figure'
            self.steps_comparison_figure(wdir, Q2, dpi)
            # self.loop_over_steps(wdir,'simple',last)
            ## 'wdir' for 'self.ratio_figure' is a dictionary that contains integer keys
            ## the keys have to be the number of the steps
            ## 'xf' of larger key values will be plotted on top of that of smaller key values
            ## 'xf' of larger key values will also be plotted with less opacity
            ## wdir = {2: './results1/step02/', 3: './results1/step03/', 4: './results1/step04/'}
            ## input file of largest key value will be loaded

        if  task == 10:
            self.msg = 'mplots.normalization_histogram'
            self.histogram_figure(wdir, Q2, dpi)
            # self.loop_over_steps(wdir,'simple',last)

    def line_figure(self, wdir, Q2, dpi):
        load_config('%s/input.py' % wdir)
        istep = core.get_istep()
        labels = load('%s/data/labels-%d.dat' % (wdir, istep))
        self.cluster = labels['cluster']
        n_replicas = len(self.cluster)
        self.best_cluster = 0

        nrows, ncols = 3, 2
        fig = py.figure(figsize = (ncols * 5.0, nrows * 3.0))
        AX = [py.subplot(nrows, ncols, cnt + 1) for cnt in range(6) if cnt != 5]

        ## PDF
        # self.xf_data = load('%s/data/pdf-%d.dat' % (wdir, istep))
        if Q2 == None: Q2 = conf['Q20']
        if Q2 == conf['Q20']:
            self.xf_data = load('%s/data/pdf-%d.dat' % (wdir, istep))
        else:
            self.xf_data = load('%s/data/pdf-%d-%f.dat' % (wdir, istep, Q2))

        print '\nplotting PDF line figure from %s' % wdir
        self.plot_lines(AX[0], 'uv', 'r', n_replicas)
        self.plot_lines(AX[0], 'dv', 'b', n_replicas)
        self.plot_lines(AX[1], 'd/u' , 'r', n_replicas)
        self.plot_lines(AX[2], 'db+ub', 'r', n_replicas)
        self.plot_lines(AX[3], 'db-ub', 'r', n_replicas)
        self.plot_lines(AX[4], 'g', 'r', n_replicas)

        AX[0].set_ylim(0.0, 0.85)
        AX[0].set_yticks([0.0, 0.2, 0.4, 0.6, 0.8])
        AX[0].set_yticklabels([r'$0$', r'$0.2$', r'$0.4$', r'$0.6$', r'$0.8$'])
        AX[0].text(0.4, 0.81, r'\boldmath$xu_v$', color = 'k', transform = AX[0].transAxes, size = 28)
        AX[0].text(0.4, 0.07, r'\boldmath$xd_v$', color = 'k', transform = AX[0].transAxes, size = 28)
        AX[1].set_ylim(0.0, 1.0)
        AX[1].set_yticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
        AX[1].set_yticklabels([r'$0$', r'', r'$0.2$', r'', r'$0.4$', r'', r'$0.6$', r'', r'$0.8$'])
        AX[1].text(0.07, 0.13, r'\boldmath$d/u$', color = 'k', transform = AX[1].transAxes, size = 28)
        AX[2].set_ylim(0.0, 0.5)
        AX[2].set_yticks([0.0, 0.1, 0.2, 0.3, 0.4])
        AX[2].set_yticklabels([r'$0$', r'', r'$0.2$', r'', r'$0.4$'])
        AX[2].text(0.07, 0.1, r'\boldmath$x(\bar{d}\!+\!\bar{u})$', color = 'k', transform = AX[2].transAxes, size = 28)
        AX[3].set_ylim(-0.05, 0.13)
        AX[3].set_yticks([-0.05, 0.0, 0.05, 0.1])
        AX[3].set_yticklabels([r'-$0.05$', r'$0$', r'$0.05$', r'$0.1$'])
        AX[3].text(0.07, 0.81, r'\boldmath$x(\bar{d}\!-\!\bar{u})$', color = 'k', transform = AX[3].transAxes, size = 28)
        AX[4].set_ylim(0.0, 2.5)
        AX[4].text(0.8, 0.8, r'\boldmath$xg$', color = 'k', transform = AX[4].transAxes, size = 31)
        AX[4].set_yticks([0.0, 0.5, 1.0, 1.5, 2.0])

        for ax in AX:
            ax.semilogx()
            ax.set_xlim(8e-3, 9e-1)
            ax.set_xticks([1e-2, 1e-1])
            ax.xaxis.set_label_coords(0.95, -0.05)
            ax.tick_params(axis = 'both', which = 'both', right = True, top = True, direction = 'in', labelsize = 20)
            ax.set_xticklabels([r'', r''])
            ax.legend(frameon = 0, loc = 'best')

        AX[0].set_ylabel(r'$xf(x)$', size = 20)
        AX[3].set_xlabel(r'\boldmath$x$', size = 30)
        AX[3].set_xticklabels([r'$0.01$', r'$0.1$'])
        AX[4].set_xlabel(r'\boldmath$x$', size = 30)
        AX[4].set_xticklabels([r'$0.01$', r'$0.1$'])

        # py.subplots_adjust(left = 0.12, bottom = 0.08, right = 0.99, top = 0.97, wspace = None, hspace = 0.2)
        py.tight_layout()
        if Q2 == conf['Q20']:
            py.savefig('%s/gallery/pdf-lines-no-strange-%d.png' % (wdir, istep), dpi = dpi)
        else:
            py.savefig('%s/gallery/pdf-lines-no-strange-%d-%f.png' % (wdir, istep, Q2), dpi = dpi)

    def groups_comparison_figure(self, wdir, Q2, dpi):
        load_config('%s/input.py' % wdir)
        istep = core.get_istep()
        labels = load('%s/data/labels-%d.dat' % (wdir, istep))
        self.cluster = labels['cluster']
        n_replicas = len(self.cluster)
        self.best_cluster = 0

        if Q2 == None: Q2 = conf['Q20']
        groups = load('%s/analysis/qpdlib/lhapdf-%f.dat' % (os.environ['FITPACK'], Q2))

        nrows, ncols = 3, 2
        fig = py.figure(figsize = (ncols * 5.0, nrows * 3.0))
        AX = [py.subplot(nrows, ncols, cnt + 1) for cnt in range(6) if cnt != 5]

        ## PDF
        if Q2 == conf['Q20']:
            self.xf_data = load('%s/data/pdf-%d.dat' % (wdir, istep))
        else:
            self.xf_data = load('%s/data/pdf-%d-%f.dat' % (wdir, istep, Q2))

        print '\nplotting PDF band figure with following groups at Q2 = %f' % Q2
        print groups['groups']

        ax = AX[0]
        self.plot_band_JAM(ax, 'uv', 'r', r'$\rm JAM$')
        for _ in groups['XF']:
            self.plot_band(ax, groups, _, 'uv', alpha = 0.2, label = r'$\rm %s$' % _)
        self.plot_band_JAM(ax, 'dv', 'r', r'$\rm JAM$')
        for _ in groups['XF']:
            self.plot_band(ax, groups, _, 'dv', alpha = 0.2, label = r'$\rm %s$' % _)
        ax.text(0.4, 0.81, r'\boldmath$xu_v$', color = 'k', transform = ax.transAxes, size = 28)
        ax.text(0.6, 0.1, r'\boldmath$xd_v$', color = 'k', transform = ax.transAxes, size = 28)
        ax.set_ylim(0.0, 0.85)
        ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8])
        ax.set_yticklabels([r'$0$', r'$0.2$', r'$0.4$', r'$0.6$', r'$0.8$'])

        ax = AX[1]
        self.plot_band_JAM(ax, 'd/u', 'r', r'$\rm JAM$')
        for _ in groups['XF']:
            self.plot_band(ax, groups, _, 'd/u', alpha = 0.2, label = r'$\rm %s$' % _)
        ax.text(0.07, 0.13, r'\boldmath$d/u$', color = 'k', transform = ax.transAxes, size = 30)
        ax.set_ylim(0.0, 0.9)
        ax.set_yticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
        ax.set_yticklabels([r'$0$', r'', r'$0.2$', r'', r'$0.4$', r'', r'$0.6$', r'', r'$0.8$'])

        ax = AX[2]
        self.plot_band_JAM(ax, 'db+ub','r', r'$\rm JAM$')
        for _ in groups['XF']:
            self.plot_band(ax, groups, _, 'db+ub', alpha = 0.2, label = r'$\rm %s$' % _)
        ax.text(0.07, 0.13, r'\boldmath$x(\bar{d}\!+\!\bar{u})$', color = 'k', transform = ax.transAxes, size = 28)
        ax.set_ylim(0.0, 0.5)
        ax.set_yticks([0.0, 0.1, 0.2, 0.3, 0.4])
        ax.set_yticklabels([r'$0$', r'', r'$0.2$', r'', r'$0.4$'])

        ax = AX[3]
        self.plot_band_JAM(ax, 'db-ub', 'r', r'$\rm JAM$')
        for _ in groups['XF']:
            self.plot_band(ax, groups, _, 'db-ub', alpha = 0.2, label = r'$\rm %s$' % _)
        ax.text(0.07, 0.8, r'\boldmath$x(\bar{d}\!-\!\bar{u})$', color = 'k', transform = ax.transAxes,size=28)
        ax.axhline(y = 0.0, c = 'k', ls = '-', lw = 1, alpha = 0.5)
        ax.set_ylim(-0.07, 0.13)
        ax.set_yticks([-0.05, 0.0, 0.05, 0.1])
        ax.set_yticklabels([r'-$0.05$', r'$0$', r'$0.05$', r'$0.1$'])

        ax = AX[4]
        self.plot_band_JAM(ax, 'g', 'r', r'$\rm JAM~Jet$')
        for _ in groups['XF']:
            if _ == 'CJ15nlo': label=r'$\rm CJ15$'
            elif _ == 'NNPDF31_nlo_as_0118': label = r'$\rm NNPDF3.1$'
            elif _ == 'MMHT2014nlo68cl': label = r'$\rm MMHT14$'
            elif _ == 'ABMP16_3_nlo': label = r'$\rm ABMP16$'
            elif _ == 'CSKK_nnlo_EIG': label = r'$\rm CSSK$'
            elif _ == 'JAM19PDF_proton_nlo': label = r'$\rm JAM19$'
            else: label = r'$\rm %s$' % _
            self.plot_band(ax, groups, _, 'g', alpha = 0.2, label = label)
        ax.text(0.8, 0.8, r'\boldmath$xg$', color = 'k', transform = ax.transAxes, size = 31)
        ax.set_ylim(0.0, 3.7)
        # ax.set_yticks([0.0, 0.5, 1.0, 1.5, 2.0])
        ax.set_yticks([0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5])
        ax.legend(frameon = 0, loc = 'best', bbox_to_anchor = [1.87, 0.9], fontsize = 15)

        for ax in AX:
            ax.semilogx()
            ax.tick_params(axis = 'both', which = 'both', right = True, top = True, direction = 'in', labelsize = 20)
            ax.xaxis.set_label_coords(0.7, -0.05)
            ax.set_xlim(8e-3, 9e-1)
            ax.set_xticks([1e-2, 1e-1])
            if ax == AX[3] or ax == AX[4]:
                ax.set_xlabel(r'\boldmath$x$', size = 30)   
                ax.set_xticklabels([r'$0.01$', r'$0.1$'])
                ax.text(0.84, -0.095, r'$0.5$', color = 'k', transform = ax.transAxes, size = 20)
            else:
                ax.set_xticklabels([r'', r''])

        py.tight_layout()
        # py.subplots_adjust(left=0.05, bottom=0.04, right=0.99, top=0.99, wspace=0.17, hspace=0.04)
        if Q2 == conf['Q20']:
            py.savefig('%s/gallery/pdf-groups-no-strange-%d.png' % (wdir, istep), dpi = dpi)
        else:
            py.savefig('%s/gallery/pdf-groups-no-strange-%d-%f.png' % (wdir, istep, Q2), dpi = dpi)

    def ratio_figure(self, wdirs, Q2, dpi):
        load_config('%s/input.py' % wdirs['current'])
        istep = core.get_istep()
        labels = load('%s/data/labels-%d.dat' % (wdirs['current'], istep))
        self.cluster = labels['cluster']
        n_replicas = len(self.cluster)
        self.best_cluster = 0
        base_step = int(re.findall(r'\d+', wdirs['base'].split('/')[2])[0])

        nrows, ncols = 4, 2
        fig = py.figure(figsize = (ncols * 5.0, nrows * 3.0))
        AX = [py.subplot(nrows, ncols, cnt+1) for cnt in range(8) if cnt != 7]

        ## PDF
        if Q2 == None: Q2 = conf['Q20']
        if Q2 == conf['Q20']:
            self.xf_data = load('%s/data/pdf-%d.dat' % (wdirs['current'], istep))
            self.base_xf_data = load('%s/data/pdf-%d.dat' % (wdirs['base'], base_step))
        else:
            self.xf_data = load('%s/data/pdf-%d-%f.dat' % (wdirs['current'], istep, Q2))
            self.base_xf_data = load('%s/data/pdf-%d-%f.dat' % (wdirs['base'], base_step, Q2))

        self.base_xf_data['mean_xf'] = {}
        for flavor in self.base_xf_data['XF']:
            self.base_xf_data['mean_xf'][flavor] = np.mean(self.base_xf_data['XF'][flavor], axis = 0)

        print '\nplotting PDF ratio figure of %s over %s at Q2 = %f' % (wdirs['current'], wdirs['base'], Q2)
        self.plot_ratios(AX[0], 'uv', 'r', n_replicas)
        self.plot_ratios(AX[0], 'dv', 'b', n_replicas)
        self.plot_ratios(AX[1], 'd/u' , 'r', n_replicas)
        self.plot_ratios(AX[2], 'db+ub', 'r', n_replicas)
        self.plot_ratios(AX[3], 'db-ub', 'r', n_replicas)
        self.plot_ratios(AX[4], 's+sb', 'r', n_replicas)
        self.plot_ratios(AX[5], 'rs', 'r', n_replicas)
        self.plot_ratios(AX[6], 'g', 'r', n_replicas)

        AX[0].text(0.4, 0.81, r'\boldmath$xu_v$', color = 'r', transform = AX[0].transAxes, size = 28)
        AX[0].text(0.4, 0.07, r'\boldmath$xd_v$', color = 'b', transform = AX[0].transAxes, size = 28)
        AX[1].text(0.07, 0.13, r'\boldmath$d/u$', color = 'k', transform = AX[1].transAxes, size = 28)
        AX[2].text(0.07, 0.1, r'\boldmath$x(\bar{d}\!+\!\bar{u})$', color = 'k', transform = AX[2].transAxes, size = 28)
        AX[3].text(0.07, 0.81, r'\boldmath$x(\bar{d}\!-\!\bar{u})$', color = 'k', transform = AX[3].transAxes, size = 28)
        AX[4].text(0.07, 0.8, r'\boldmath$x(s\!+\!\bar{s})$', color = 'k', transform = AX[4].transAxes, size = 28)
        AX[5].text(0.05, 0.78, r'\boldmath$R_s$', color = 'k', transform = AX[5].transAxes, size = 30)
        AX[6].text(0.07, 0.8, r'\boldmath$xg$', color = 'k', transform = AX[6].transAxes, size = 31)

        for ax in AX:
            ax.semilogx()
            ax.set_xlim(8e-3, 9e-1)
            ax.set_ylim(0.8, 1.2)
            ax.set_yticks([0.9, 1.0, 1.1])
            ax.set_yticklabels([r'$0.9$', r'$1.0$', r'$1.1$'])
            ax.set_xticks([1e-2, 1e-1])
            ax.xaxis.set_label_coords(0.95, -0.05)
            ax.tick_params(axis = 'both', which = 'both', right = True, top = True, direction = 'in', labelsize = 20)
            ax.set_xticklabels([r'', r''])
            ax.axhline(1.0, color = 'k', linestyle = ':')
            ax.legend(frameon = 0, loc = 'best')

        AX[0].set_ylabel(r'$xf(x)$', size = 20)
        AX[5].set_xlabel(r'\boldmath$x$', size = 30)
        AX[5].set_xticklabels([r'$0.01$', r'$0.1$'])
        AX[6].set_xlabel(r'\boldmath$x$', size = 30)
        AX[6].set_xticklabels([r'$0.01$', r'$0.1$'])

        # py.subplots_adjust(left = 0.12, bottom = 0.08, right = 0.99, top = 0.97, wspace = None, hspace = 0.2)
        py.tight_layout()
        if Q2 == conf['Q20']:
            py.savefig('%s/gallery/pdf-ratio-no-strange-%d.png' % (wdirs['current'], istep), dpi = dpi)
        else:
            py.savefig('%s/gallery/pdf-ratio-no-strange-%d-%f.png' % (wdirs['current'], istep, Q2), dpi = dpi)

    def steps_comparison_figure(self, wdirs, Q2, dpi):
        step_keys = sorted(wdirs.keys())
        load_config('%s/input.py' % wdirs[step_keys[-1]])
        istep = core.get_istep()

        nrows, ncols = 3, 2
        fig = py.figure(figsize = (ncols * 5.0, nrows * 3.0))
        AX = [py.subplot(nrows, ncols, cnt + 1) for cnt in range(6) if cnt != 5]
        colors = ['coral', 'grey', 'wheat', 'blueviolet', 'brown', 'olivedrab', 'deepskyblue']
        # colors = ['k', 'r', 'g', 'Yellow', 'b']

        ## PDF
        self.xf_data = {}
        self.cluster = {}
        self.best_cluster = {}
        labels = {}
        if Q2 == None: Q2 = conf['Q20']
        if Q2 == conf['Q20']:
            for step_key in step_keys:
                self.xf_data[step_key] = load('%s/data/pdf-%d.dat' % (wdirs[step_key], step_key))
                labels[step_key] = load('%s/data/labels-%d.dat' % (wdirs[step_key], step_key))
                self.cluster[step_key] = labels[step_key]['cluster']
                self.best_cluster[step_key] = 0
        else:
            for step_key in step_keys:
                self.xf_data[step_key] = load('%s/data/pdf-%d-%f.dat' % (wdirs[step_key], step_key, Q2))
                labels[step_key] = load('%s/data/labels-%d.dat' % (wdirs[step_key], step_key))
                self.cluster[step_key] = labels[step_key]['cluster']
                self.best_cluster[step_key] = 0

        print '\nplotting PDF band figure with following steps at Q2 = %f' % Q2
        print [wdirs[_] for _ in wdirs]

        ax = AX[0]
        for i in range(len(step_keys)):
            step = step_keys[i]
            # alpha = (1.0 / len(step_keys)) * (len(step_keys) - i)
            alpha = 0.5
            color = colors[i]
            self.plot_band_steps(ax, step, 'uv', color, alpha)
            self.plot_band_steps(ax, step, 'dv', color, alpha)
        ax.text(0.4, 0.81, r'\boldmath$xu_v$', color = 'k', transform = ax.transAxes, size = 28)
        ax.text(0.6, 0.1, r'\boldmath$xd_v$', color = 'k', transform = ax.transAxes, size = 28)
        ax.set_ylim(0.0, 0.85)
        ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8])
        ax.set_yticklabels([r'$0$', r'$0.2$', r'$0.4$', r'$0.6$', r'$0.8$'])

        ax = AX[1]
        for i in range(len(step_keys)):
            step = step_keys[i]
            # alpha = (1.0 / len(step_keys)) * (len(step_keys) - i)
            alpha = 0.5
            color = colors[i]
            self.plot_band_steps(ax, step, 'd/u', color, alpha)
        ax.text(0.07, 0.13, r'\boldmath$d/u$', color = 'k', transform = ax.transAxes, size = 30)
        ax.set_ylim(0.0, 0.9)
        ax.set_yticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
        ax.set_yticklabels([r'$0$', r'', r'$0.2$', r'', r'$0.4$', r'', r'$0.6$', r'', r'$0.8$'])

        ax = AX[2]
        for i in range(len(step_keys)):
            step = step_keys[i]
            # alpha = (1.0 / len(step_keys)) * (len(step_keys) - i)
            alpha = 0.5
            color = colors[i]
            self.plot_band_steps(ax, step, 'db+ub', color, alpha)
        ax.text(0.07, 0.13, r'\boldmath$x(\bar{d}\!+\!\bar{u})$', color = 'k', transform = ax.transAxes, size = 28)
        ax.set_ylim(0.0, 0.5)
        ax.set_yticks([0.0, 0.1, 0.2, 0.3, 0.4])
        ax.set_yticklabels([r'$0$', r'', r'$0.2$', r'', r'$0.4$'])

        ax = AX[3]
        for i in range(len(step_keys)):
            step = step_keys[i]
            # alpha = (1.0 / len(step_keys)) * (len(step_keys) - i)
            alpha = 0.5
            color = colors[i]
            self.plot_band_steps(ax, step, 'db-ub', color, alpha)
        ax.text(0.07, 0.8, r'\boldmath$x(\bar{d}\!-\!\bar{u})$', color = 'k', transform = ax.transAxes,size=28)
        ax.axhline(y = 0.0, c = 'k', ls = '-', lw = 1, alpha = 0.5)
        ax.set_ylim(-0.07, 0.13)
        ax.set_yticks([-0.05, 0.0, 0.05, 0.1])
        ax.set_yticklabels([r'-$0.05$', r'$0$', r'$0.05$', r'$0.1$'])

        ax = AX[4]
        for i in range(len(step_keys)):
            step = step_keys[i]
            # alpha = (1.0 / len(step_keys)) * (len(step_keys) - i)
            alpha = 0.5
            color = colors[i]
            label = r'$\mathrm{step~%02d}$' % step
            self.plot_band_steps(ax, step, 'g', color, alpha, label)
        ax.text(0.8, 0.8, r'\boldmath$xg$', color = 'k', transform = ax.transAxes, size = 31)
        ax.set_ylim(0.0, 2.5)
        ax.set_yticks([0.0, 0.5, 1.0, 1.5, 2.0])
        ax.legend(frameon = 0, loc = 'best', bbox_to_anchor = [1.87, 0.9], fontsize = 20)

        for ax in AX:
            ax.semilogx()
            ax.tick_params(axis = 'both', which = 'both', right = True, top = True, direction = 'in', labelsize = 20)
            ax.xaxis.set_label_coords(0.7, -0.05)
            ax.set_xlim(8e-3, 9e-1)
            ax.set_xticks([1e-2, 1e-1])
            if ax == AX[3] or ax == AX[4]:
                ax.set_xlabel(r'\boldmath$x$', size = 30)   
                ax.set_xticklabels([r'$0.01$', r'$0.1$'])
                ax.text(0.84, -0.095, r'$0.5$', color = 'k', transform = ax.transAxes, size = 20)
            else:
                ax.set_xticklabels([r'', r''])

        py.tight_layout()
        # py.subplots_adjust(left=0.05, bottom=0.04, right=0.99, top=0.99, wspace=0.17, hspace=0.04)
        if Q2 == conf['Q20']:
            py.savefig('%s/gallery/pdf-steps-no-strange-%d.png' % (wdirs[step_keys[-1]], istep), dpi = dpi)
        else:
            py.savefig('%s/gallery/pdf-steps-no-strange-%d-%f.png' % (wdirs[step_keys[-1]], istep, Q2), dpi = dpi)

    def histogram_figure(self, wdir, Q2, dpi):
        ## plot histogram for normalizaton of all flavors except sea
        ## because first and second shape sea are fixed to be the same
        load_config('%s/input.py' % wdir)
        istep = core.get_istep()
        # jar = load('%s/data/jar-%d.dat' % (wdir, istep))
        labels = load('%s/data/labels-%d.dat' % (wdir, istep))
        self.cluster = labels['cluster']
        n_replicas = len(self.cluster)
        self.best_cluster = 0

        if Q2 == None: Q2 = conf['Q20']
        if Q2 == conf['Q20']:
            self.norm_data = load('%s/data/pdf-normalization-%d.dat' % (wdir, istep))
        else:
            self.norm_data= load('%s/data/pdf-normalization-%d-%f.dat' % (wdir, istep, Q2))

        # replicas = core.get_replicas(wdir)
        # core.mod_conf(istep, replicas[0]) ## set conf as specified in istep

        # if 'pdf' not in conf['steps'][istep]['active distributions']:
        #     print('pdf-proton not in active distribution')
        #     return

        # resman = RESMAN(nworkers = 1, parallel = False, datasets = False)
        # parman = resman.parman
        # parman.order = replicas[0]['order'][istep] ## make sure 'parman' uses the same order as all the replicas do
        # # print parman.order

        # pdf = conf['pdf']

        # core.mod_conf(istep, replicas[i_replica - 1])
        # parman.set_new_params(replicas[i_replica - 1]['params'][istep], initial = True)
        # for name, value in pdf.params.iteritems():
        #     print '%7s: %.5e, %.5e, %.5e, %.5e, %.5e' % (name, value[0], value[1], value[2], value[3], value[4])

        nrows, ncols = 4, 2
        fig = py.figure(figsize = (ncols * 5.0, nrows * 3.0))
        # AX = [py.subplot(nrows, ncols, cnt) for cnt in range(1, 9)]
        AX = [py.subplot(nrows, ncols, cnt+1) for cnt in range(8) if cnt != 7]

        ## PDF
        print '\nplotting PDF normalization histogram from %s' % wdir
        self.n_histogram(AX[0], 'uv', n_replicas, r'\boldmath$u_v$')
        self.n_histogram(AX[1], 'dv', n_replicas, r'\boldmath$d_v$')
        self.n_histogram(AX[2], 'ub' , n_replicas, r'\boldmath$u_b$')
        self.n_histogram(AX[3], 'db', n_replicas, r'\boldmath$d_b$')
        self.n_histogram(AX[4], 's', n_replicas, r'\boldmath$s$')
        self.n_histogram(AX[5], 'sb', n_replicas, r'\boldmath$s_b$')
        self.n_histogram(AX[6], 'g', n_replicas, r'\boldmath$g$')
        # self.n_histogram(AX[7], 'sea', n_replicas, r'\boldmath$sea$')

        for ax in AX:
            # ax.semilogx()
            # ax.set_xlim(8e-3, 9e-1)
            # ax.set_xticks([1e-2, 1e-1])
            # ax.set_xlabel(r'\boldmath$x$',size=15)
            # ax.xaxis.set_label_coords(0.95, -0.05)
            # ax.tick_params(axis = 'both', which = 'both', right = True, top = True, direction = 'in', labelsize = 20)
            # ax.set_xticklabels([r'', r''])
            ax.legend(frameon = 0, loc = 'best', fontsize = 17)

        ## log scale x axis has to be disabled for the time being
        ## because we need to try negative values for normalization of second shape
        # AX[2].semilogx()
        # AX[3].semilogx()
        # AX[4].semilogx()
        # AX[5].semilogx()

        # py.subplots_adjust(left = 0.12, bottom = 0.08, right = 0.99, top = 0.97, wspace = None, hspace = 0.2)
        py.tight_layout()
        if Q2 == conf['Q20']:
            py.savefig('%s/gallery/pdf-normalization-histogram-%d.png' % (wdir, istep), dpi = dpi)
        else:
            py.savefig('%s/gallery/pdf-normalization-histogram-%d-%f.png' % (wdir, istep, Q2), dpi = dpi)

class JET(PDF_PLOT_CORE):

    def __init__(self, task, wdir, Q2 = None, dpi = 200):

        if  task == 1:
            self.msg = 'mplots.line_figure'
            self.line_figure(wdir, Q2, dpi)
            # self.loop_over_steps(wdir,'simple',last)

        if  task == 2:
            self.msg = 'mplots.groups_comparison_figure'
            self.groups_comparison_figure(wdir, Q2, dpi)
            # self.loop_over_steps(wdir,'simple',last)

        if  task == 3:
            self.msg = 'mplots.ratio_figure'
            self.ratio_figure(wdir, Q2, dpi)
            # self.loop_over_steps(wdir,'simple',last)
            ## 'wdir' for 'self.ratio_figure' is a dictionary that contains keys 'current' and 'base'
            ## 'self.ratio_figure' will plot ratio of 'current' 'xf' values over 'base' 'xf' values
            ## wdir = {'base': './results1/step04/', 'current': './results1/step06/'}
            ## input file of current step will be loaded

        if  task == 4:
            self.msg = 'mplots.steps_comparison_figure'
            self.steps_comparison_figure(wdir, Q2, dpi)
            # self.loop_over_steps(wdir,'simple',last)
            ## 'wdir' for 'self.ratio_figure' is a dictionary that contains integer keys
            ## the keys have to be the number of the steps
            ## 'xf' of larger key values will be plotted on top of that of smaller key values
            ## 'xf' of larger key values will also be plotted with less opacity
            ## wdir = {2: './results1/step02/', 3: './results1/step03/', 4: './results1/step04/'}
            ## input file of largest key value will be loaded

        if  task == 10:
            self.msg = 'mplots.normalization_histogram'
            self.histogram_figure(wdir, Q2, dpi)
            # self.loop_over_steps(wdir,'simple',last)

    def line_figure(self, wdir, Q2, dpi):
        load_config('%s/input.py' % wdir)
        istep = core.get_istep()
        labels = load('%s/data/labels-%d.dat' % (wdir, istep))
        self.cluster = labels['cluster']
        n_replicas = len(self.cluster)
        self.best_cluster = 0

        nrows, ncols = 3, 2
        fig = py.figure(figsize = (ncols * 5.0, nrows * 3.0))
        AX = [py.subplot(nrows, ncols, cnt + 1) for cnt in range(6) if cnt != 5]

        ## PDF
        # self.xf_data = load('%s/data/pdf-%d.dat' % (wdir, istep))
        if Q2 == None: Q2 = conf['Q20']
        if Q2 == conf['Q20']:
            self.xf_data = load('%s/data/pdf-%d.dat' % (wdir, istep))
        else:
            self.xf_data = load('%s/data/pdf-%d-%f.dat' % (wdir, istep, Q2))

        print '\nplotting PDF line figure from %s' % wdir
        self.plot_lines(AX[0], 'uv', 'r', n_replicas)
        self.plot_lines(AX[0], 'dv', 'b', n_replicas)
        self.plot_lines(AX[1], 'd/u' , 'r', n_replicas)
        self.plot_lines(AX[2], 'db+ub', 'r', n_replicas)
        self.plot_lines(AX[3], 'db-ub', 'r', n_replicas)
        self.plot_lines(AX[4], 'g', 'r', n_replicas)

        AX[0].set_ylim(0.0, 0.85)
        AX[0].set_yticks([0.0, 0.2, 0.4, 0.6, 0.8])
        AX[0].set_yticklabels([r'$0$', r'$0.2$', r'$0.4$', r'$0.6$', r'$0.8$'])
        AX[0].text(0.4, 0.81, r'\boldmath$xu_v$', color = 'k', transform = AX[0].transAxes, size = 28)
        AX[0].text(0.4, 0.07, r'\boldmath$xd_v$', color = 'k', transform = AX[0].transAxes, size = 28)
        AX[1].set_ylim(0.0, 1.0)
        AX[1].set_yticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
        AX[1].set_yticklabels([r'$0$', r'', r'$0.2$', r'', r'$0.4$', r'', r'$0.6$', r'', r'$0.8$'])
        AX[1].text(0.07, 0.13, r'\boldmath$d/u$', color = 'k', transform = AX[1].transAxes, size = 28)
        AX[2].set_ylim(0.0, 0.5)
        AX[2].set_yticks([0.0, 0.1, 0.2, 0.3, 0.4])
        AX[2].set_yticklabels([r'$0$', r'', r'$0.2$', r'', r'$0.4$'])
        AX[2].text(0.07, 0.1, r'\boldmath$x(\bar{d}\!+\!\bar{u})$', color = 'k', transform = AX[2].transAxes, size = 28)
        AX[3].set_ylim(-0.05, 0.13)
        AX[3].set_yticks([-0.05, 0.0, 0.05, 0.1])
        AX[3].set_yticklabels([r'-$0.05$', r'$0$', r'$0.05$', r'$0.1$'])
        AX[3].text(0.07, 0.81, r'\boldmath$x(\bar{d}\!-\!\bar{u})$', color = 'k', transform = AX[3].transAxes, size = 28)
        AX[4].set_ylim(0.0, 2.5)
        AX[4].text(0.8, 0.8, r'\boldmath$xg$', color = 'k', transform = AX[4].transAxes, size = 31)
        AX[4].set_yticks([0.0, 0.5, 1.0, 1.5, 2.0])

        for ax in AX:
            ax.semilogx()
            ax.set_xlim(8e-3, 9e-1)
            ax.set_xticks([1e-2, 1e-1])
            ax.xaxis.set_label_coords(0.95, -0.05)
            ax.tick_params(axis = 'both', which = 'both', right = True, top = True, direction = 'in', labelsize = 20)
            ax.set_xticklabels([r'', r''])
            ax.legend(frameon = 0, loc = 'best')

        AX[0].set_ylabel(r'$xf(x)$', size = 20)
        AX[3].set_xlabel(r'\boldmath$x$', size = 30)
        AX[3].set_xticklabels([r'$0.01$', r'$0.1$'])
        AX[4].set_xlabel(r'\boldmath$x$', size = 30)
        AX[4].set_xticklabels([r'$0.01$', r'$0.1$'])

        # py.subplots_adjust(left = 0.12, bottom = 0.08, right = 0.99, top = 0.97, wspace = None, hspace = 0.2)
        py.tight_layout()
        if Q2 == conf['Q20']:
            py.savefig('%s/gallery/pdf-lines-jet-%d.png' % (wdir, istep), dpi = dpi)
        else:
            py.savefig('%s/gallery/pdf-lines-jet-%d-%f.png' % (wdir, istep, Q2), dpi = dpi)

    def groups_comparison_figure(self, wdir, Q2, dpi):
        load_config('%s/input.py' % wdir)
        istep = core.get_istep()
        labels = load('%s/data/labels-%d.dat' % (wdir, istep))
        self.cluster = labels['cluster']
        n_replicas = len(self.cluster)
        self.best_cluster = 0

        if Q2 == None: Q2 = conf['Q20']
        groups = load('%s/analysis/qpdlib/lhapdf-%f.dat' % (os.environ['FITPACK'], Q2))

        nrows, ncols = 3, 2
        fig = py.figure(figsize = (ncols * 5.0, nrows * 3.0))
        AX = [py.subplot(nrows, ncols, cnt + 1) for cnt in range(6) if cnt != 5]

        ## PDF
        if Q2 == conf['Q20']:
            self.xf_data = load('%s/data/pdf-%d.dat' % (wdir, istep))
        else:
            self.xf_data = load('%s/data/pdf-%d-%f.dat' % (wdir, istep, Q2))

        print '\nplotting PDF band figure with following groups at Q2 = %f' % Q2
        print groups['groups']

        jet_distinction = {'alpha': 0.4, 'zorder': 1}

        ax = AX[0]
        self.plot_band_JAM(ax, 'uv', 'r', r'$\rm JAM~Jet$', distinction = jet_distinction)
        for _ in groups['XF']:
            self.plot_band(ax, groups, _, 'uv', alpha = 0.2, label = r'$\rm %s$' % _)
        self.plot_band_JAM(ax, 'dv', 'r', r'$\rm JAM~Jet$', distinction = jet_distinction)
        for _ in groups['XF']:
            self.plot_band(ax, groups, _, 'dv', alpha = 0.2, label = r'$\rm %s$' % _)
        ax.text(0.4, 0.81, r'\boldmath$xu_v$', color = 'k', transform = ax.transAxes, size = 28)
        ax.text(0.6, 0.1, r'\boldmath$xd_v$', color = 'k', transform = ax.transAxes, size = 28)
        ax.set_ylim(0.0, 0.85)
        ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8])
        ax.set_yticklabels([r'$0$', r'$0.2$', r'$0.4$', r'$0.6$', r'$0.8$'])

        ax = AX[1]
        self.plot_band_JAM(ax, 'd/u', 'r', r'$\rm JAM~Jet$', distinction = jet_distinction)
        for _ in groups['XF']:
            self.plot_band(ax, groups, _, 'd/u', alpha = 0.2, label = r'$\rm %s$' % _)
        ax.text(0.07, 0.13, r'\boldmath$d/u$', color = 'k', transform = ax.transAxes, size = 30)
        ax.set_ylim(0.0, 0.9)
        ax.set_yticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
        ax.set_yticklabels([r'$0$', r'', r'$0.2$', r'', r'$0.4$', r'', r'$0.6$', r'', r'$0.8$'])

        ax = AX[2]
        self.plot_band_JAM(ax, 'db+ub','r', r'$\rm JAM~Jet$', distinction = jet_distinction)
        for _ in groups['XF']:
            self.plot_band(ax, groups, _, 'db+ub', alpha = 0.2, label = r'$\rm %s$' % _)
        ax.text(0.07, 0.13, r'\boldmath$x(\bar{d}\!+\!\bar{u})$', color = 'k', transform = ax.transAxes, size = 28)
        ax.set_ylim(0.0, 0.5)
        ax.set_yticks([0.0, 0.1, 0.2, 0.3, 0.4])
        ax.set_yticklabels([r'$0$', r'', r'$0.2$', r'', r'$0.4$'])

        ax = AX[3]
        self.plot_band_JAM(ax, 'db-ub', 'r', r'$\rm JAM~Jet$', distinction = jet_distinction)
        for _ in groups['XF']:
            self.plot_band(ax, groups, _, 'db-ub', alpha = 0.2, label = r'$\rm %s$' % _)
        ax.text(0.07, 0.8, r'\boldmath$x(\bar{d}\!-\!\bar{u})$', color = 'k', transform = ax.transAxes,size=28)
        ax.axhline(y = 0.0, c = 'k', ls = '-', lw = 1, alpha = 0.5)
        ax.set_ylim(-0.07, 0.13)
        ax.set_yticks([-0.05, 0.0, 0.05, 0.1])
        ax.set_yticklabels([r'-$0.05$', r'$0$', r'$0.05$', r'$0.1$'])

        ax = AX[4]
        self.plot_band_JAM(ax, 'g', 'r', r'$\rm JAM~Jet$', distinction = jet_distinction)
        for _ in groups['XF']:
            if _ == 'CJ15nlo': label=r'$\rm CJ15$'
            elif _ == 'NNPDF31_nlo_as_0118': label = r'$\rm NNPDF3.1$'
            elif _ == 'MMHT2014nlo68cl': label = r'$\rm MMHT14$'
            elif _ == 'ABMP16_3_nlo': label = r'$\rm ABMP16$'
            elif _ == 'CSKK_nnlo_EIG': label = r'$\rm CSSK$'
            elif _ == 'JAM19PDF_proton_nlo': label = r'$\rm JAM19$'
            else: label = r'$\rm %s$' % _
            self.plot_band(ax, groups, _, 'g', alpha = 0.2, label = label)
        ax.text(0.8, 0.8, r'\boldmath$xg$', color = 'k', transform = ax.transAxes, size = 31)
        ax.set_ylim(0.0, 3.7)
        # ax.set_yticks([0.0, 0.5, 1.0, 1.5, 2.0])
        ax.set_yticks([0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5])
        ax.legend(frameon = 0, loc = 'best', bbox_to_anchor = [1.87, 0.9], fontsize = 15)

        for ax in AX:
            ax.semilogx()
            ax.tick_params(axis = 'both', which = 'both', right = True, top = True, direction = 'in', labelsize = 20)
            ax.xaxis.set_label_coords(0.7, -0.05)
            ax.set_xlim(8e-3, 9e-1)
            ax.set_xticks([1e-2, 1e-1])
            if ax == AX[3] or ax == AX[4]:
                ax.set_xlabel(r'\boldmath$x$', size = 30)   
                ax.set_xticklabels([r'$0.01$', r'$0.1$'])
                ax.text(0.84, -0.095, r'$0.5$', color = 'k', transform = ax.transAxes, size = 20)
            else:
                ax.set_xticklabels([r'', r''])

        # fig.suptitle(r'$Q^2 ~=~ %.2f$' % Q2)
        py.tight_layout()
        py.subplots_adjust(left=0.05, bottom=0.07, right=0.99, top=0.99, wspace=0.17, hspace=0.2)
        if Q2 == conf['Q20']:
            py.savefig('%s/gallery/pdf-groups-jet-%d.png' % (wdir, istep), dpi = dpi)
        else:
            py.savefig('%s/gallery/pdf-groups-jet-%d-%f.png' % (wdir, istep, Q2), dpi = dpi)

    def ratio_figure(self, wdirs, Q2, dpi):
        load_config('%s/input.py' % wdirs['current'])
        istep = core.get_istep()
        labels = load('%s/data/labels-%d.dat' % (wdirs['current'], istep))
        self.cluster = labels['cluster']
        n_replicas = len(self.cluster)
        self.best_cluster = 0
        base_step = int(re.findall(r'\d+', wdirs['base'].split('/')[2])[0])

        nrows, ncols = 4, 2
        fig = py.figure(figsize = (ncols * 5.0, nrows * 3.0))
        AX = [py.subplot(nrows, ncols, cnt+1) for cnt in range(8) if cnt != 7]

        ## PDF
        if Q2 == None: Q2 = conf['Q20']
        if Q2 == conf['Q20']:
            self.xf_data = load('%s/data/pdf-%d.dat' % (wdirs['current'], istep))
            self.base_xf_data = load('%s/data/pdf-%d.dat' % (wdirs['base'], base_step))
        else:
            self.xf_data = load('%s/data/pdf-%d-%f.dat' % (wdirs['current'], istep, Q2))
            self.base_xf_data = load('%s/data/pdf-%d-%f.dat' % (wdirs['base'], base_step, Q2))

        self.base_xf_data['mean_xf'] = {}
        for flavor in self.base_xf_data['XF']:
            self.base_xf_data['mean_xf'][flavor] = np.mean(self.base_xf_data['XF'][flavor], axis = 0)

        print '\nplotting PDF ratio figure of %s over %s at Q2 = %f' % (wdirs['current'], wdirs['base'], Q2)
        self.plot_ratios(AX[0], 'uv', 'r', n_replicas)
        self.plot_ratios(AX[0], 'dv', 'b', n_replicas)
        self.plot_ratios(AX[1], 'd/u' , 'r', n_replicas)
        self.plot_ratios(AX[2], 'db+ub', 'r', n_replicas)
        self.plot_ratios(AX[3], 'db-ub', 'r', n_replicas)
        self.plot_ratios(AX[4], 's+sb', 'r', n_replicas)
        self.plot_ratios(AX[5], 'rs', 'r', n_replicas)
        self.plot_ratios(AX[6], 'g', 'r', n_replicas)

        AX[0].text(0.4, 0.81, r'\boldmath$xu_v$', color = 'r', transform = AX[0].transAxes, size = 28)
        AX[0].text(0.4, 0.07, r'\boldmath$xd_v$', color = 'b', transform = AX[0].transAxes, size = 28)
        AX[1].text(0.07, 0.13, r'\boldmath$d/u$', color = 'k', transform = AX[1].transAxes, size = 28)
        AX[2].text(0.07, 0.1, r'\boldmath$x(\bar{d}\!+\!\bar{u})$', color = 'k', transform = AX[2].transAxes, size = 28)
        AX[3].text(0.07, 0.81, r'\boldmath$x(\bar{d}\!-\!\bar{u})$', color = 'k', transform = AX[3].transAxes, size = 28)
        AX[4].text(0.07, 0.8, r'\boldmath$x(s\!+\!\bar{s})$', color = 'k', transform = AX[4].transAxes, size = 28)
        AX[5].text(0.05, 0.78, r'\boldmath$R_s$', color = 'k', transform = AX[5].transAxes, size = 30)
        AX[6].text(0.07, 0.8, r'\boldmath$xg$', color = 'k', transform = AX[6].transAxes, size = 31)

        for ax in AX:
            ax.semilogx()
            ax.set_xlim(8e-3, 9e-1)
            ax.set_ylim(0.8, 1.2)
            ax.set_yticks([0.9, 1.0, 1.1])
            ax.set_yticklabels([r'$0.9$', r'$1.0$', r'$1.1$'])
            ax.set_xticks([1e-2, 1e-1])
            ax.xaxis.set_label_coords(0.95, -0.05)
            ax.tick_params(axis = 'both', which = 'both', right = True, top = True, direction = 'in', labelsize = 20)
            ax.set_xticklabels([r'', r''])
            ax.axhline(1.0, color = 'k', linestyle = ':')
            ax.legend(frameon = 0, loc = 'best')

        AX[0].set_ylabel(r'$xf(x)$', size = 20)
        AX[5].set_xlabel(r'\boldmath$x$', size = 30)
        AX[5].set_xticklabels([r'$0.01$', r'$0.1$'])
        AX[6].set_xlabel(r'\boldmath$x$', size = 30)
        AX[6].set_xticklabels([r'$0.01$', r'$0.1$'])

        # py.subplots_adjust(left = 0.12, bottom = 0.08, right = 0.99, top = 0.97, wspace = None, hspace = 0.2)
        py.tight_layout()
        if Q2 == conf['Q20']:
            py.savefig('%s/gallery/pdf-ratio-jet-%d.png' % (wdirs['current'], istep), dpi = dpi)
        else:
            py.savefig('%s/gallery/pdf-ratio-jet-%d-%f.png' % (wdirs['current'], istep, Q2), dpi = dpi)

    def steps_comparison_figure(self, wdirs, Q2, dpi):
        step_keys = sorted(wdirs.keys())
        load_config('%s/input.py' % wdirs[step_keys[-1]])
        istep = core.get_istep()

        nrows, ncols = 1, 1
        fig = py.figure(figsize = (ncols * 5.0, nrows * 3.0))
        # AX = [py.subplot(nrows, ncols, cnt + 1) for cnt in range(6) if cnt != 5]
        ax = py.subplot(nrows, ncols, 1)
        colors = ['coral', 'grey', 'wheat', 'blueviolet', 'brown', 'olivedrab', 'deepskyblue']
        # colors = ['k', 'r', 'g', 'Yellow', 'b']

        ## PDF
        self.xf_data = {}
        self.cluster = {}
        self.best_cluster = {}
        labels = {}
        if Q2 == None: Q2 = conf['Q20']
        if Q2 == conf['Q20']:
            for step_key in step_keys:
                self.xf_data[step_key] = load('%s/data/pdf-%d.dat' % (wdirs[step_key], step_key))
                labels[step_key] = load('%s/data/labels-%d.dat' % (wdirs[step_key], step_key))
                self.cluster[step_key] = labels[step_key]['cluster']
                self.best_cluster[step_key] = 0
        else:
            for step_key in step_keys:
                self.xf_data[step_key] = load('%s/data/pdf-%d-%f.dat' % (wdirs[step_key], step_key, Q2))
                labels[step_key] = load('%s/data/labels-%d.dat' % (wdirs[step_key], step_key))
                self.cluster[step_key] = labels[step_key]['cluster']
                self.best_cluster[step_key] = 0

        print '\nplotting PDF band figure with following steps at Q2 = %f' % Q2
        print [wdirs[_] for _ in wdirs]

        for i in range(len(step_keys)):
            step = step_keys[i]
            # alpha = (1.0 / len(step_keys)) * (len(step_keys) - i)
            alpha = 0.5
            color = colors[i]
            label = r'$\mathrm{step~%02d}$' % step
            self.plot_band_steps(ax, step, 'g', color, alpha, label)
        ax.text(0.3, 0.3, r'\boldmath$xg$', color = 'k', transform = ax.transAxes, size = 31)
        ax.set_ylim(0.0, 2.5)
        ax.set_yticks([0.0, 0.5, 1.0, 1.5, 2.0])
        ax.legend(frameon = 0, loc = 'upper right', fontsize = 10)

        ax.semilogx()
        ax.tick_params(axis = 'both', which = 'both', right = True, top = True, direction = 'in', labelsize = 20)
        ax.xaxis.set_label_coords(0.7, -0.05)
        ax.set_xlim(8e-3, 9e-1)
        ax.set_xticks([1e-2, 1e-1])
        ax.set_xlabel(r'\boldmath$x$', size = 30)   
        ax.set_xticklabels([r'$0.01$', r'$0.1$'])
        ax.text(0.84, -0.095, r'$0.5$', color = 'k', transform = ax.transAxes, size = 20)

        py.tight_layout()
        # py.subplots_adjust(left=0.05, bottom=0.04, right=0.99, top=0.99, wspace=0.17, hspace=0.04)
        if Q2 == conf['Q20']:
            py.savefig('%s/gallery/pdf-steps-jet-%d.png' % (wdirs[step_keys[-1]], istep), dpi = dpi)
        else:
            py.savefig('%s/gallery/pdf-steps-jet-%d-%f.png' % (wdirs[step_keys[-1]], istep, Q2), dpi = dpi)

    def histogram_figure(self, wdir, Q2, dpi):
        ## plot histogram for normalizaton of all flavors except sea
        ## because first and second shape sea are fixed to be the same
        load_config('%s/input.py' % wdir)
        istep = core.get_istep()
        # jar = load('%s/data/jar-%d.dat' % (wdir, istep))
        labels = load('%s/data/labels-%d.dat' % (wdir, istep))
        self.cluster = labels['cluster']
        n_replicas = len(self.cluster)
        self.best_cluster = 0

        if Q2 == None: Q2 = conf['Q20']
        if Q2 == conf['Q20']:
            self.norm_data = load('%s/data/pdf-normalization-%d.dat' % (wdir, istep))
        else:
            self.norm_data= load('%s/data/pdf-normalization-%d-%f.dat' % (wdir, istep, Q2))

        # replicas = core.get_replicas(wdir)
        # core.mod_conf(istep, replicas[0]) ## set conf as specified in istep

        # if 'pdf' not in conf['steps'][istep]['active distributions']:
        #     print('pdf-proton not in active distribution')
        #     return

        # resman = RESMAN(nworkers = 1, parallel = False, datasets = False)
        # parman = resman.parman
        # parman.order = replicas[0]['order'][istep] ## make sure 'parman' uses the same order as all the replicas do
        # # print parman.order

        # pdf = conf['pdf']

        # core.mod_conf(istep, replicas[i_replica - 1])
        # parman.set_new_params(replicas[i_replica - 1]['params'][istep], initial = True)
        # for name, value in pdf.params.iteritems():
        #     print '%7s: %.5e, %.5e, %.5e, %.5e, %.5e' % (name, value[0], value[1], value[2], value[3], value[4])

        nrows, ncols = 4, 2
        fig = py.figure(figsize = (ncols * 5.0, nrows * 3.0))
        # AX = [py.subplot(nrows, ncols, cnt) for cnt in range(1, 9)]
        AX = [py.subplot(nrows, ncols, cnt+1) for cnt in range(8) if cnt != 7]

        ## PDF
        print '\nplotting PDF normalization histogram from %s' % wdir
        self.n_histogram(AX[0], 'uv', n_replicas, r'\boldmath$u_v$')
        self.n_histogram(AX[1], 'dv', n_replicas, r'\boldmath$d_v$')
        self.n_histogram(AX[2], 'ub' , n_replicas, r'\boldmath$u_b$')
        self.n_histogram(AX[3], 'db', n_replicas, r'\boldmath$d_b$')
        self.n_histogram(AX[4], 's', n_replicas, r'\boldmath$s$')
        self.n_histogram(AX[5], 'sb', n_replicas, r'\boldmath$s_b$')
        self.n_histogram(AX[6], 'g', n_replicas, r'\boldmath$g$')
        # self.n_histogram(AX[7], 'sea', n_replicas, r'\boldmath$sea$')

        for ax in AX:
            # ax.semilogx()
            # ax.set_xlim(8e-3, 9e-1)
            # ax.set_xticks([1e-2, 1e-1])
            # ax.set_xlabel(r'\boldmath$x$',size=15)
            # ax.xaxis.set_label_coords(0.95, -0.05)
            # ax.tick_params(axis = 'both', which = 'both', right = True, top = True, direction = 'in', labelsize = 20)
            # ax.set_xticklabels([r'', r''])
            ax.legend(frameon = 0, loc = 'best', fontsize = 17)

        ## log scale x axis has to be disabled for the time being
        ## because we need to try negative values for normalization of second shape
        # AX[2].semilogx()
        # AX[3].semilogx()
        # AX[4].semilogx()
        # AX[5].semilogx()

        # py.subplots_adjust(left = 0.12, bottom = 0.08, right = 0.99, top = 0.97, wspace = None, hspace = 0.2)
        py.tight_layout()
        if Q2 == conf['Q20']:
            py.savefig('%s/gallery/pdf-normalization-histogram-%d.png' % (wdir, istep), dpi = dpi)
        else:
            py.savefig('%s/gallery/pdf-normalization-histogram-%d-%f.png' % (wdir, istep, Q2), dpi = dpi)

if __name__ == '__main__':
    pass
