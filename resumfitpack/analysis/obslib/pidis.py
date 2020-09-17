#!/usr/bin/env python
import sys, os
import numpy as np
import copy
from subprocess import Popen, PIPE, STDOUT
import pandas as pd
from scipy import stats
from scipy.stats import gaussian_kde
from scipy.integrate import quad
from mpmath import fp
import tempfile

## matplotlib
import matplotlib
matplotlib.use('Agg')
import pylab as py
import matplotlib.gridspec as gridspec
from  matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Times-Roman']})
rc('text',usetex=True)
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
from matplotlib.ticker import FixedLocator

## from fitpack tools
from tools.tools     import load, save, checkdir, lprint, tex
from tools.config    import conf, load_config

## from fitpack fitlib
from fitlib.resman import RESMAN

## from fitpack analysis
from analysis.corelib import core
from analysis.corelib import classifier

def generate_maps():
    maps = {}

    maps['experiment'] = {}
    maps['experiment']['EMC']=tex('EMC')
    maps['experiment']['SMC']=tex('SMC')
    maps['experiment']['SLACE80E130']=tex('E80/E130')
    maps['experiment']['SLAC(E142)']=tex('E142')
    maps['experiment']['SLAC(E143)']=tex('E143')
    maps['experiment']['SLAC(E154)']=tex('E154')
    maps['experiment']['SLAC(E155)']=tex('E155')
    maps['experiment']['SLAC(E155x)']=tex('E155x')
    maps['experiment']['HERMES']=tex('HERMES')
    maps['experiment']['COMPASS']=tex('COMPASS')
    maps['experiment']['JLabHA(E06014)']=tex('E06014')
    maps['experiment']['JLabHA(E99117)']=tex('E99117')
    maps['experiment']['JLabHB(EG1DVCS)']=tex('eg1')+'-'+tex('DVCS')
    maps['experiment']['JLabHB(EG1b)']=tex('eg1b')

    maps['target'] = {}
    maps['target']['proton']='p'
    maps['target']['neutron']='n'
    maps['target']['helium']='He'
    maps['target']['deuteron']='d'

    maps['observable'] = {}
    maps['observable']['A1']  = r'$A_1^{%s}$'
    maps['observable']['A2']  = r'$A_2^{%s}$'
    maps['observable']['Apa'] = r'$A_\parallel^{%s}$'
    maps['observable']['Ape'] = r'$A_\perp^{%s}$'
    maps['observable']['Atpe']= r'$\tilde{A}_\perp^{%s}$'

    maps['color'] = {}
    maps['color']['p']='r'
    maps['color']['He']='b'
    maps['color']['d']='g'

    return maps

def pick_subsets(data):
    set_headers = {}
    for dataset in data:
        set_heater_temp = [k for k in data[dataset].keys() if 'set' in k]
        if any(set_heater_temp): set_headers[dataset] = set_heater_temp[0]

    subset_bins = {}
    for dataset in data:
        if dataset in set_headers:
            subset_bins[dataset] = {}
            set_header = set_headers[dataset]
            subsets = [int(i) for i in data[dataset][set_header]] ## some sets can be absent due to 'Q2' cut and 'W2' cut
            subsets = sorted(set(subsets))
            for subset in subsets:
                subset_temp = []
                for i in range(len(data[dataset][set_header])):
                    if data[dataset][set_header][i] == subset: subset_temp.append(i)
                subset_bins[dataset][subset] = subset_temp
        else:
            subset_bins[dataset] = {}
            subset_bins[dataset][0] = range(len(data[dataset]['X']))

    return subset_bins

def data_Apa_proton(wdir, data, indices, subset_bins, log_list, kc, istep, dpi):
    n_columns, n_rows = 4, 5
    py.figure(figsize = (n_columns * 2.5, n_rows * 1.3))
    ax_s = {}
    count = 1
    for index in indices:
        count += 1
        ax = py.subplot(n_rows, n_columns, count)
        correlation = [_ for _ in data[index] if '_c' in _]
        i_sets = np.sort(subset_bins[index].keys())
        for i_set in i_sets:
            idx = subset_bins[index][i_set]
            if (data[index]['col'] == 'SLACE80E130') and (i_set == 5): continue
            if len(idx) <= 2: continue
        print any(correlation)
    return

def data_Ape_proton(wdir, data, indices, subset_bins, log_list, kc, istep, dpi):
    print indices
    for index in indices:
        correlation = [_ for _ in data[index] if '_c' in _]
        print any(correlation)
    return

def data_Apa_proton_DVCS(wdir, data, indices, subset_bins, log_list, kc, istep, dpi):
    print indices
    for index in indices:
        correlation = [_ for _ in data[index] if '_c' in _]
        print any(correlation)
    return

def data_Apa_proton_EG1b(wdir, data, indices, subset_bins, log_list, kc, istep, dpi):
    print indices
    for index in indices:
        correlation = [_ for _ in data[index] if '_c' in _]
        print any(correlation)
    return

def data_Apa_deuteron(wdir, data, indices, subset_bins, log_list, kc, istep, dpi):
    print indices
    for index in indices:
        correlation = [_ for _ in data[index] if '_c' in _]
        print any(correlation)
    return

def data_Ape_deuteron(wdir, data, indices, subset_bins, log_list, kc, istep, dpi):
    print indices
    for index in indices:
        correlation = [_ for _ in data[index] if '_c' in _]
        print any(correlation)
    return

def data_Apa_deuteron_DVCS(wdir, data, indices, subset_bins, log_list, kc, istep, dpi):
    print indices
    for index in indices:
        correlation = [_ for _ in data[index] if '_c' in _]
        print any(correlation)
    return

def data_Apa_deuteron_EG1b(wdir, data, indices, subset_bins, log_list, kc, istep, dpi):
    print indices
    for index in indices:
        correlation = [_ for _ in data[index] if '_c' in _]
        print any(correlation)
    return

def data_helium(wdir, data, indices, subset_bins, log_list, kc, istep, dpi):
    print indices
    for index in indices:
        correlation = [_ for _ in data[index] if '_c' in _]
        print any(correlation)
    return

def make_figure(wdir, task, kc, dpi = 200):

    if task == 1:
        print('\nplotting PIDIS data from %s' % (wdir))
    elif task == 2:
        print('\nplotting PIDIS data over theory from %s' % (wdir))

    load_config('%s/input.py' % wdir)
    istep = core.get_istep()
    replicas = core.get_replicas(wdir)
    core.mod_conf(istep, replicas[0]) #--set conf as specified in istep

    predictions = load('%s/data/predictions-%d.dat' % (wdir, istep))
    if 'pidis' not in predictions['reactions']:
        print('PIDIS is not in data file')
        return
    labels  = load('%s/data/labels-%d.dat' % (wdir, istep))
    cluster = labels['cluster']

    data = predictions['reactions']['pidis']

    for idx in data:
        predictions = copy.copy(data[idx]['prediction-rep'])
        del data[idx]['prediction-rep']
        del data[idx]['residuals-rep']
        for ic in range(kc.nc[istep]):
            predictions_ic = [predictions[i] for i in range(len(predictions)) if cluster[i] == ic]
            data[idx]['thy-%d' % ic] = np.mean(predictions_ic, axis = 0)
            data[idx]['dthy-%d' % ic] = np.std(predictions_ic, axis = 0) ** 0.5

    maps = generate_maps()
    subset_bins = pick_subsets(data)
    log_list = []

    if task == 1:
        indices = {}
        for idx in data:
            observable = data[idx]['obs'][0]
            target = data[idx]['target'][0]
            collaboration = data[idx]['col'][0]
            if (observable == 'A1') or (observable == 'Apa'):
                if target == 'p':
                    if collaboration.startswith('JLab') != True:
                        if 'Apa_proton' not in indices:
                            indices['Apa_proton'] = []
                        indices['Apa_proton'].append(idx)
                    elif collaboration.startswith('JLabHB(EG1DVCS)') == True:
                        if 'Apa_proton_DVCS' not in indices:
                            indices['Apa_proton_DVCS'] = []
                        indices['Apa_proton_DVCS'].append(idx)
                    elif collaboration.startswith('JLabHB(EG1b)') == True:
                        if 'Apa_proton_EG1b' not in indices:
                            indices['Apa_proton_EG1b'] = []
                        indices['Apa_proton_EG1b'].append(idx)
                elif target == 'd':
                    if collaboration.startswith('JLab') != True:
                        if 'Apa_deuteron' not in indices:
                            indices['Apa_deuteron'] = []
                        indices['Apa_deuteron'].append(idx)
                    elif collaboration.startswith('JLabHB(EG1DVCS)') == True:
                        if 'Apa_deuteron_DVCS' not in indices:
                            indices['Apa_deuteron_DVCS'] = []
                        indices['Apa_deuteron_DVCS'].append(idx)
                    elif collaboration.startswith('JLabHB(EG1b)') == True:
                        if 'Apa_deuteron_EG1b' not in indices:
                            indices['Apa_deuteron_EG1b'] = []
                        indices['Apa_deuteron_EG1b'].append(idx)
            elif (observable == 'A2') or (observable == 'Ape') or (observable == 'Atpe'):
                if target == 'p':
                    if collaboration.startswith('JLab') != True:
                        if 'Ape_proton' not in indices:
                            indices['Ape_proton'] = []
                        indices['Ape_proton'].append(idx)
                elif target == 'd':
                    if collaboration.startswith('JLab') != True:
                        if 'Ape_deuteron' not in indices:
                            indices['Ape_deuteron'] = []
                        indices['Ape_deuteron'].append(idx)

            if target == 'h':
                if 'helium' not in indices:
                    indices['helium'] = []
                indices['helium'].append(idx)

        for key in indices:
            if key == 'Apa_proton':
                data_Apa_proton(wdir, data, indices[key], subset_bins, log_list, kc, istep, dpi)
            elif key == 'Apa_proton_DVCS':
                data_Apa_proton_DVCS(wdir, data, indices[key], subset_bins, log_list, kc, istep, dpi)
            elif key == 'Apa_proton_EG1b':
                data_Apa_proton_EG1b(wdir, data, indices[key], subset_bins, log_list, kc, istep, dpi)
            elif key == 'Apa_deuteron':
                data_Apa_deuteron(wdir, data, indices[key], subset_bins, log_list, kc, istep, dpi)
            elif key == 'Apa_deuteron_DVCS':
                data_Apa_deuteron_DVCS(wdir, data, indices[key], subset_bins, log_list, kc, istep, dpi)
            elif key == 'Apa_deuteron_EG1b':
                data_Apa_deuteron_EG1b(wdir, data, indices[key], subset_bins, log_list, kc, istep, dpi)
            elif key == 'Ape_proton':
                data_Ape_proton(wdir, data, indices[key], subset_bins, log_list, kc, istep, dpi)
            elif key == 'Ape_deuteron':
                data_Ape_deuteron(wdir, data, indices[key], subset_bins, log_list, kc, istep, dpi)
            elif key == 'helium':
                data_helium(wdir, data, indices[key], subset_bins, log_list, kc, istep, dpi)

    elif task == 2:
        plot_data_on_theory(wdir, data, kc, istep, dpi)

    return

if __name__ == "__main__":
    pass
