#!/usr/bin/env python
import os,sys
import subprocess
import numpy as np
import scipy as sp
import pandas as pd
import copy

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
# from pdfcalc  import PDFCALC

class PPDF_PLOT_CORE:


    def plot_lines(self, ax, flav, c, n_replicas, lab = ''): ## plot all the PDF replicas as separate lines
        line_width = 10.0 / float(n_replicas)
        if line_width > 1.0: line_width = 1.0
        X = self.xf_data['X']
        n_replicas = len(self.xf_data['XF'][flav])
        for i in range(n_replicas):
            if self.cluster[i] != self.best_cluster: continue
            ax.plot(X, self.xf_data['XF'][flav][i], color = c, label = lab, linewidth = line_width)

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
        if flav == 'up':
            ax.hist(self.norm_data['first'][flav + '1'], n_replicas, histtype = 'bar', color = 'r', label = label.rsplit('$', 1)[0] + '^{(1)}$')
            ax.hist(self.norm_data['second'][flav + '2'], n_replicas, histtype = 'bar', color = 'b', label = label.rsplit('$', 1)[0] + '^{(2)}$')
        else:
            ax.hist(self.norm_data['first'][flav + '1'], n_replicas, histtype = 'bar', color = 'r', label = label)
            ax.hist(self.norm_data['second'][flav + '2'], n_replicas, histtype = 'bar', color = 'b')

class PLOTS(PPDF_PLOT_CORE):

    def __init__(self, task, wdir, Q2 = None, dpi = 200):

        if  task == 1:
            self.msg = 'mplots.line_figure'
            self.line_figure(wdir, Q2, dpi)
            # self.loop_over_steps(wdir,'simple',last)

        if  task == 4:
            self.msg = 'mplots.steps_comparison_figure'
            self.steps_comparison_figure(wdir, Q2, dpi)
            # self.loop_over_steps(wdir,'simple',last)
            ## 'wdir' for 'self.ratio_figure' is a dictionary that contains integer keys
            ## the keys have to be the number of the steps
            ## 'xf' of larger key values will be plotted on top of that of smaller key values
            ## 'xf' of larger key values will also be plotted with less opacity
            ## wdir = {6: './results1/step06/', 7: './results1/step07/', 9: './results1/step09/'}
            ## input file of largest key value will be loaded

        if  task == 10:
            self.msg = 'mplots.normalization_histogram'
            self.histogram_figure(wdir, Q2, dpi)
            # self.loop_over_steps(wdir,'simple',last)

    def line_figure(self, wdir, Q2, dpi):
        ## this block may contain something not useful
        load_config('%s/input.py' % wdir)
        istep = core.get_istep()
        # jar = load('%s/data/jar-%d.dat' % (wdir, istep))
        labels = load('%s/data/labels-%d.dat' % (wdir, istep))
        self.cluster = labels['cluster']
        n_replicas = len(self.cluster)
        self.best_cluster = 0
        ## this block may contain something not useful

        nrows, ncols = 2, 2
        fig = py.figure(figsize = (ncols * 5.0, nrows * 3.0))
        AX = [py.subplot(nrows, ncols, cnt + 1) for cnt in range(4)]

        ## polarized PDF
        if Q2 == None: Q2 = conf['Q20']
        if Q2 == conf['Q20']:
            self.xf_data = load('%s/data/ppdf-%d.dat' % (wdir, istep))
        else:
            self.xf_data = load('%s/data/ppdf-%d-%f.dat' % (wdir, istep, Q2))

        print '\nplotting polarized PDF line figure from %s at Q2 = %f' % (wdir, Q2)
        self.plot_lines(AX[0], 'up', 'r', n_replicas)
        self.plot_lines(AX[1], 'dp' , 'r', n_replicas)
        self.plot_lines(AX[2], 'sp', 'r', n_replicas)
        self.plot_lines(AX[3], 'g', 'r', n_replicas)

        AX[0].set_ylim(0.0, 0.5)
        AX[0].set_yticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5])
        AX[0].set_yticklabels([r'$0.1$', r'$0.2$', r'$0.3$', r'$0.4$', r'$0.5$'])
        AX[0].text(0.6, 0.3, r'\boldmath$x\Delta u^+$', color = 'k', transform = AX[0].transAxes, size = 28)
        AX[1].set_ylim(-0.2, 0.0)
        AX[1].set_yticks([-0.05, -0.1, -0.15])
        AX[1].set_yticklabels([r'$-0.05$', r'$-0.10$', r'$-0.15$'])
        AX[1].text(0.05, 0.15, r'\boldmath$x\Delta d^+$', color = 'k', transform = AX[1].transAxes, size = 28)
        AX[2].set_ylim(-0.1, 0.05)
        AX[2].set_yticks([-0.08, -0.04, 0.0, 0.04])
        AX[2].set_yticklabels([r'$-0.08$', r'$-0.04$', r'$0.00$', r'$0.04$'])
        AX[2].text(0.1, 0.7, r'\boldmath$x\Delta s^+$', color = 'k', transform = AX[2].transAxes, size = 28)
        AX[3].set_ylim(-0.11, 0.3)
        AX[3].set_yticks([-0.1, 0.0, 0.1, 0.2])
        AX[3].set_yticklabels([r'$-0.1$', r'$0.0$', r'$0.1$', r'$0.2$'])
        AX[3].text(0.7, 0.8, r'\boldmath$x\Delta g$', color = 'k', transform = AX[3].transAxes, size = 28)

        for ax in AX:
            ax.semilogx()
            ax.set_xlim(8e-3, 9e-1)
            ax.set_xticks([1e-2, 1e-1])
            # ax.xaxis.set_label_coords(0.95, -0.05)
            ax.tick_params(axis = 'both', which = 'both', right = True, top = True, direction = 'in', labelsize = 20)
            ax.set_xticklabels([r'', r''])
            ax.legend(frameon = 0, loc = 'best')

        AX[2].axhline(0.0, color = 'k', linestyle = ':')
        AX[3].axhline(0.0, color = 'k', linestyle = ':')
        # AX[0].set_ylabel(r'$xf(x)$', size = 20)
        # AX[5].set_xlabel(r'\boldmath$x$', size = 30)
        # AX[5].set_xticklabels([r'$0.01$', r'$0.1$'])
        # AX[6].set_xlabel(r'\boldmath$x$', size = 30)
        AX[2].set_xticks([0.01, 0.1, 0.5, 0.8])
        AX[2].set_xticklabels([r'$0.01$', r'$0.1$', r'$0.5$', r'$0.8$'])
        AX[3].set_xticks([0.01, 0.1, 0.5, 0.8])
        AX[3].set_xticklabels([r'$0.01$', r'$0.1$', r'$0.5$', r'$0.8$'])

        py.subplots_adjust(left = 0.12, bottom = 0.08, right = 0.99, top = 0.97, wspace = None, hspace = 0.2)
        py.tight_layout()
        if Q2 == conf['Q20']:
            py.savefig('%s/gallery/ppdf-lines-%d.png' % (wdir, istep), dpi = dpi)
        else:
            py.savefig('%s/gallery/ppdf-lines-%d-%f.png' % (wdir, istep, Q2), dpi = dpi)

    def steps_comparison_figure(self, wdirs, Q2, dpi):
        step_keys = sorted(wdirs.keys())
        load_config('%s/input.py' % wdirs[step_keys[-1]])
        istep = core.get_istep()

        nrows, ncols = 2, 2
        fig = py.figure(figsize = (ncols * 5.0, nrows * 3.0))
        AX = [py.subplot(nrows, ncols, cnt + 1) for cnt in range(4)]
        colors = ['k', 'r', 'g', 'Yellow', 'b']

        ## PPDF
        self.xf_data = {}
        self.cluster = {}
        self.best_cluster = {}
        labels = {}
        if Q2 == None: Q2 = conf['Q20']
        if Q2 == conf['Q20']:
            for step_key in step_keys:
                self.xf_data[step_key] = load('%s/data/ppdf-%d.dat' % (wdirs[step_key], step_key))
                labels[step_key] = load('%s/data/labels-%d.dat' % (wdirs[step_key], step_key))
                self.cluster[step_key] = labels[step_key]['cluster']
                self.best_cluster[step_key] = 0
        else:
            for step_key in step_keys:
                self.xf_data[step_key] = load('%s/data/ppdf-%d-%f.dat' % (wdirs[step_key], step_key, Q2))
                labels[step_key] = load('%s/data/labels-%d.dat' % (wdirs[step_key], step_key))
                self.cluster[step_key] = labels[step_key]['cluster']
                self.best_cluster[step_key] = 0

        print '\nplotting PPDF band figure with following steps at Q2 = %f' % Q2
        print [wdirs[_] for _ in wdirs]

        ax = AX[0]
        for i in range(len(step_keys)):
            step = step_keys[i]
            # alpha = (1.0 / len(step_keys)) * (len(step_keys) - i)
            alpha = 0.5
            color = colors[i]
            self.plot_band_steps(ax, step, 'uv', color, alpha)
            self.plot_band_steps(ax, step, 'dv', color, alpha)
        ax.set_ylim(0.0, 0.5)
        ax.set_yticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5])
        ax.set_yticklabels([r'$0.1$', r'$0.2$', r'$0.3$', r'$0.4$', r'$0.5$'])
        ax.text(0.6, 0.3, r'\boldmath$x\Delta u^+$', color = 'k', transform = AX[0].transAxes, size = 28)

        ax = AX[1]
        for i in range(len(step_keys)):
            step = step_keys[i]
            # alpha = (1.0 / len(step_keys)) * (len(step_keys) - i)
            alpha = 0.5
            color = colors[i]
            self.plot_band_steps(ax, step, 'd/u', color, alpha)
        ax.set_ylim(-0.2, 0.0)
        ax.set_yticks([-0.05, -0.1, -0.15])
        ax.set_yticklabels([r'$-0.05$', r'$-0.10$', r'$-0.15$'])
        ax.text(0.05, 0.15, r'\boldmath$x\Delta d^+$', color = 'k', transform = AX[1].transAxes, size = 28)

        ax = AX[2]
        for i in range(len(step_keys)):
            step = step_keys[i]
            # alpha = (1.0 / len(step_keys)) * (len(step_keys) - i)
            alpha = 0.5
            color = colors[i]
            self.plot_band_steps(ax, step, 'db+ub', color, alpha)
        ax.set_ylim(-0.1, 0.05)
        ax.set_yticks([-0.08, -0.04, 0.0, 0.04])
        ax.set_yticklabels([r'$-0.08$', r'$-0.04$', r'$0.00$', r'$0.04$'])
        ax.text(0.1, 0.7, r'\boldmath$x\Delta s^+$', color = 'k', transform = AX[2].transAxes, size = 28)

        ax = AX[3]
        for i in range(len(step_keys)):
            step = step_keys[i]
            # alpha = (1.0 / len(step_keys)) * (len(step_keys) - i)
            alpha = 0.5
            color = colors[i]
            self.plot_band_steps(ax, step, 'db-ub', color, alpha)
        ax.set_ylim(-0.11, 0.3)
        ax.set_yticks([-0.1, 0.0, 0.1, 0.2])
        ax.set_yticklabels([r'$-0.1$', r'$0.0$', r'$0.1$', r'$0.2$'])
        ax.text(0.7, 0.8, r'\boldmath$x\Delta g$', color = 'k', transform = AX[3].transAxes, size = 28)

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
            ax.set_xlim(8e-3, 9e-1)
            ax.set_xticks([1e-2, 1e-1])
            # ax.xaxis.set_label_coords(0.95, -0.05)
            ax.tick_params(axis = 'both', which = 'both', right = True, top = True, direction = 'in', labelsize = 20)
            ax.set_xticklabels([r'', r''])
            ax.xaxis.set_label_coords(0.7, -0.05)
            ax.legend(frameon = 0, loc = 'best')

        AX[2].axhline(0.0, color = 'k', linestyle = ':')
        AX[3].axhline(0.0, color = 'k', linestyle = ':')
        # AX[0].set_ylabel(r'$xf(x)$', size = 20)
        # AX[5].set_xlabel(r'\boldmath$x$', size = 30)
        # AX[5].set_xticklabels([r'$0.01$', r'$0.1$'])
        # AX[6].set_xlabel(r'\boldmath$x$', size = 30)
        AX[2].set_xticks([0.01, 0.1, 0.5, 0.8])
        AX[2].set_xticklabels([r'$0.01$', r'$0.1$', r'$0.5$', r'$0.8$'])
        AX[3].set_xticks([0.01, 0.1, 0.5, 0.8])
        AX[3].set_xticklabels([r'$0.01$', r'$0.1$', r'$0.5$', r'$0.8$'])

        py.tight_layout()
        py.subplots_adjust(left=0.05, bottom=0.04, right=0.99, top=0.99, wspace=0.17, hspace=0.04)
        if Q2 == conf['Q20']:
            py.savefig('%s/gallery/ppdf-steps-%d.png' % (wdirs[step_keys[-1]], istep), dpi = dpi)
        else:
            py.savefig('%s/gallery/ppdf-steps-%d-%f.png' % (wdirs[step_keys[-1]], istep, Q2), dpi = dpi)

    def histogram_figure(self, wdir, Q2, dpi = 200):
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
            self.norm_data = load('%s/data/ppdf-normalization-%d.dat' % (wdir, istep))
        else:
            self.norm_data= load('%s/data/ppdf-normalization-%d-%f.dat' % (wdir, istep, Q2))

        # replicas = core.get_replicas(wdir)
        # core.mod_conf(istep, replicas[0]) ## set conf as specified in istep

        # if 'ppdf' not in conf['steps'][istep]['active distributions']:
        #     print('ppdf-proton not in active distribution')
        #     return

        # resman = RESMAN(nworkers = 1, parallel = False, datasets = False)
        # parman = resman.parman
        # parman.order = replicas[0]['order'][istep] ## make sure 'parman' uses the same order as all the replicas do
        # # print parman.order

        # ppdf = conf['ppdf']

        # core.mod_conf(istep, replicas[i_replica - 1])
        # parman.set_new_params(replicas[i_replica - 1]['params'][istep], initial = True)
        # for name, value in ppdf.params.iteritems():
        #     print '%7s: %.5e, %.5e, %.5e, %.5e, %.5e' % (name, value[0], value[1], value[2], value[3], value[4])

        nrows, ncols = 4, 2
        fig = py.figure(figsize = (ncols * 5.0, nrows * 3.0))
        # AX = [py.subplot(nrows, ncols, cnt) for cnt in range(1, 9)]
        AX = [py.subplot(nrows, ncols, cnt+1) for cnt in range(8) if cnt != 7]

        ## PDF
        print '\nplotting PDF normalization histogram from %s at Q2 = %f' % (wdir, Q2)
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
            py.savefig('%s/gallery/ppdf-normalization-histogram-%d.png' % (wdir, istep), dpi = dpi)
        else:
            py.savefig('%s/gallery/ppdf-normalization-histogram-%d-%f.png' % (wdir, istep, Q2), dpi = dpi)

class PJET(PPDF_PLOT_CORE):

    def __init__(self, task, wdir, Q2 = None, dpi = 200):

        if  task == 1:
            self.msg = 'mplots.line_figure'
            self.line_figure(wdir, Q2, dpi)
            # self.loop_over_steps(wdir,'simple',last)

        if  task == 4:
            self.msg = 'mplots.steps_comparison_figure'
            self.steps_comparison_figure(wdir, Q2, dpi)
            # self.loop_over_steps(wdir,'simple',last)
            ## 'wdir' for 'self.ratio_figure' is a dictionary that contains integer keys
            ## the keys have to be the number of the steps
            ## 'xf' of larger key values will be plotted on top of that of smaller key values
            ## 'xf' of larger key values will also be plotted with less opacity
            ## wdir = {6: './results1/step06/', 7: './results1/step07/', 9: './results1/step09/'}
            ## input file of largest key value will be loaded

        if  task == 10:
            self.msg = 'mplots.normalization_histogram'
            self.histogram_figure(wdir, Q2, dpi)
            # self.loop_over_steps(wdir,'simple',last)

    def line_figure(self, wdir, Q2, dpi):
        ## this block may contain something not useful
        load_config('%s/input.py' % wdir)
        istep = core.get_istep()
        # jar = load('%s/data/jar-%d.dat' % (wdir, istep))
        labels = load('%s/data/labels-%d.dat' % (wdir, istep))
        self.cluster = labels['cluster']
        n_replicas = len(self.cluster)
        self.best_cluster = 0
        ## this block may contain something not useful

        nrows, ncols = 1, 1
        fig = py.figure(figsize = (ncols * 5.0, nrows * 3.0))
        ax = py.subplot(nrows, ncols, 1)

        ## polarized PDF
        if Q2 == None: Q2 = conf['Q20']
        if Q2 == conf['Q20']:
            self.xf_data = load('%s/data/ppdf-%d.dat' % (wdir, istep))
        else:
            self.xf_data = load('%s/data/ppdf-%d-%f.dat' % (wdir, istep, Q2))

        print '\nplotting polarized PDF line figure from %s at Q2 = %f' % (wdir, Q2)
        self.plot_lines(ax, 'g', 'r', n_replicas)

        ax.set_ylim(-0.11, 0.3)
        ax.set_yticks([-0.1, 0.0, 0.1, 0.2])
        ax.set_yticklabels([r'$-0.1$', r'$0.0$', r'$0.1$', r'$0.2$'])
        ax.text(0.7, 0.8, r'\boldmath$x\Delta g$', color = 'k', transform = AX[3].transAxes, size = 28)

        ax.semilogx()
        ax.set_xlim(8e-3, 9e-1)
        ax.set_xticks([1e-2, 1e-1])
        ax.tick_params(axis = 'both', which = 'both', right = True, top = True, direction = 'in', labelsize = 20)
        ax.set_xticklabels([r'', r''])
        ax.legend(frameon = 0, loc = 'best')

        ax.axhline(0.0, color = 'k', linestyle = ':')
        ax.set_xticks([0.01, 0.1, 0.5, 0.8])
        ax.set_xticklabels([r'$0.01$', r'$0.1$', r'$0.5$', r'$0.8$'])

        # py.subplots_adjust(left = 0.12, bottom = 0.08, right = 0.99, top = 0.97, wspace = None, hspace = 0.2)
        py.tight_layout()
        if Q2 == conf['Q20']:
            py.savefig('%s/gallery/ppdf-lines-pjet-%d.png' % (wdir, istep), dpi = dpi)
        else:
            py.savefig('%s/gallery/ppdf-lines-pjet-%d-%f.png' % (wdir, istep, Q2), dpi = dpi)

    def steps_comparison_figure(self, wdirs, Q2, dpi):
        step_keys = sorted(wdirs.keys())
        load_config('%s/input.py' % wdirs[step_keys[-1]])
        istep = core.get_istep()

        nrows, ncols = 1, 1
        fig = py.figure(figsize = (ncols * 5.0, nrows * 3.0))
        ax = py.subplot(nrows, ncols, 1)
        colors = ['y', 'r', 'b', 'g', 'm']

        ## PPDF
        self.xf_data = {}
        self.cluster = {}
        self.best_cluster = {}
        labels = {}
        if Q2 == None: Q2 = conf['Q20']
        if Q2 == conf['Q20']:
            for step_key in step_keys:
                self.xf_data[step_key] = load('%s/data/ppdf-%d.dat' % (wdirs[step_key], step_key))
                labels[step_key] = load('%s/data/labels-%d.dat' % (wdirs[step_key], step_key))
                self.cluster[step_key] = labels[step_key]['cluster']
                self.best_cluster[step_key] = 0
        else:
            for step_key in step_keys:
                self.xf_data[step_key] = load('%s/data/ppdf-%d-%f.dat' % (wdirs[step_key], step_key, Q2))
                labels[step_key] = load('%s/data/labels-%d.dat' % (wdirs[step_key], step_key))
                self.cluster[step_key] = labels[step_key]['cluster']
                self.best_cluster[step_key] = 0

        print '\nplotting PPDF band figure with following steps at Q2 = %f' % Q2
        print [wdirs[_] for _ in wdirs]

        for i in range(len(step_keys)):
            step = step_keys[i]
            # alpha = (1.0 / len(step_keys)) * (len(step_keys) - i)
            alpha = 0.5
            color = colors[i]
            label = r'$\mathrm{step~%02d}$' % step
            self.plot_band_steps(ax, step, 'g', color, alpha, label)
        ax.set_ylim(-0.11, 0.3)
        ax.set_yticks([-0.1, 0.0, 0.1, 0.2])
        ax.set_yticklabels([r'$-0.1$', r'$0.0$', r'$0.1$', r'$0.2$'])
        # ax.text(0.1, 0.1, r'\boldmath$x\Delta g$', color = 'k', transform = ax.transAxes, size = 28)
        ax.legend(frameon = 0, loc = 'lower left', fontsize = 7)

        ax.semilogx()
        ax.set_xlim(8e-3, 9e-1)
        ax.set_xticks([1e-2, 1e-1])
        ax.tick_params(axis = 'both', which = 'both', right = True, top = True, direction = 'in', labelsize = 20)
        ax.set_xticklabels([r'', r''])

        ax.axhline(0.0, color = 'k', linestyle = ':')
        ax.set_xticks([0.01, 0.1, 0.5, 0.8])
        ax.set_xticklabels([r'$0.01$', r'$0.1$', r'$0.5$', r'$0.8$'])
        ax.set_ylabel(r'$x \Delta g$', size = 15)
        ax.set_xlabel(r'\boldmath$x$', size = 15)
        ax.xaxis.set_label_coords(1.05, 0.0)

        py.tight_layout()
        # py.subplots_adjust(left=0.05, bottom=0.04, right=0.99, top=0.99, wspace=0.17, hspace=0.04)
        if Q2 == conf['Q20']:
            py.savefig('%s/gallery/ppdf-steps-pjet-%d.png' % (wdirs[step_keys[-1]], istep), dpi = dpi)
        else:
            py.savefig('%s/gallery/ppdf-steps-pjet-%d-%f.png' % (wdirs[step_keys[-1]], istep, Q2), dpi = dpi)

    def histogram_figure(self, wdir, Q2, dpi = 200):
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
            self.norm_data = load('%s/data/ppdf-normalization-%d.dat' % (wdir, istep))
        else:
            self.norm_data= load('%s/data/ppdf-normalization-%d-%f.dat' % (wdir, istep, Q2))

        # replicas = core.get_replicas(wdir)
        # core.mod_conf(istep, replicas[0]) ## set conf as specified in istep

        # if 'ppdf' not in conf['steps'][istep]['active distributions']:
        #     print('ppdf-proton not in active distribution')
        #     return

        # resman = RESMAN(nworkers = 1, parallel = False, datasets = False)
        # parman = resman.parman
        # parman.order = replicas[0]['order'][istep] ## make sure 'parman' uses the same order as all the replicas do
        # # print parman.order

        # ppdf = conf['ppdf']

        # core.mod_conf(istep, replicas[i_replica - 1])
        # parman.set_new_params(replicas[i_replica - 1]['params'][istep], initial = True)
        # for name, value in ppdf.params.iteritems():
        #     print '%7s: %.5e, %.5e, %.5e, %.5e, %.5e' % (name, value[0], value[1], value[2], value[3], value[4])

        nrows, ncols = 4, 2
        fig = py.figure(figsize = (ncols * 5.0, nrows * 3.0))
        # AX = [py.subplot(nrows, ncols, cnt) for cnt in range(1, 9)]
        AX = [py.subplot(nrows, ncols, cnt+1) for cnt in range(8) if cnt != 7]

        ## PDF
        print '\nplotting PDF normalization histogram from %s at Q2 = %f' % (wdir, Q2)
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
            py.savefig('%s/gallery/ppdf-normalization-histogram-%d.png' % (wdir, istep), dpi = dpi)
        else:
            py.savefig('%s/gallery/ppdf-normalization-histogram-%d-%f.png' % (wdir, istep, Q2), dpi = dpi)

if __name__ == '__main__':
    pass
