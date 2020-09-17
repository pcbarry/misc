import os
import numpy as np

import matplotlib
matplotlib.use('Agg')
import pylab as py


#--from tools
from tools.tools     import checkdir,save,load
import tools.config
from tools.config    import load_config, conf, options
from tools.inputmod  import INPUTMOD

#--from local
from analysis.corelib import core
from analysis.corelib import classifier


def plot_params(wdir,dist,kc,histogram=False):
    load_config('%s/input.py'%wdir)
    istep=core.get_istep()
    replicas = core.get_replicas(wdir)
    core.mod_conf(istep,replicas[0])

    labels = load('%s/data/labels-%d.dat' % (wdir, istep))
    #clusters = labels['cluster']
    #colors = labels['cluster_colors']
    clusters,colors,nc,cluster_order = classifier.get_clusters(wdir,istep,kc)

    _orders = replicas[0]['order'][istep]

    #--get correct orders from dist
    orders  = []
    indices    = []
    for i in range(len(_orders)):
        if _orders[i][1] != dist: continue
        orders.append(_orders[i][2])
        indices.append(i)
    if len(orders) == 0:
        print('%s is a passive distribution' % dist)
        return
    orders = np.array(orders)
    indices = np.array(indices)
    ## this part is necessary for sorting individual parameters
    sorting_indices = np.argsort(orders)
    orders = orders[sorting_indices]
    indices = indices[sorting_indices]

    flavors = sorted(list(set([_.split(' ')[0] for _ in orders])))
    print('\nmaking scatter plot for %s parameters of the following flavors' % dist)
    print(flavors)
    n_parameter_max = 0
    for flavor in flavors:
        count = 0
        for order in orders:
            if order.split(' ')[0] == flavor:
                count += 1
        if count > n_parameter_max:
            n_parameter_max = count

    #--get correct params from dist
    params = np.zeros((len(orders),len(replicas)))
    for i in range(len(orders)):
        for j in range(len(replicas)):
            params[i][j] = replicas[j]['params'][istep][indices[i]]

    #--create plot with enough space for # of parameters
    # nrows, ncols = np.ceil(len(orders)/ float(n_parameter_max)), n_parameter_max
    nrows, ncols = len(flavors), n_parameter_max
    fig = py.figure(figsize = (ncols * 7, nrows * 4))
    X = np.linspace(1, len(replicas), len(replicas))

    #--create plot
    i_row = 0
    i_ax = 0
    for flavor in flavors:
        i_ax = i_row * n_parameter_max
        for i in range(len(orders)):
            if orders[i].split(' ')[0] == flavor:
                ax = py.subplot(nrows, ncols, i_ax + 1)
                ## if following chunck causes errors, comment below out and use instead
                ## ax.set_title('%s'%(orders[i]), size=20)
                if 'uv' in orders[i].split(' ')[0]:
                    shape = int(orders[i].split(' ')[0][-1])
                    title = r'$u_v^{\left( %d \right)}~%s$' % (shape, orders[i].split(' ')[1])
                elif 'ubv' in orders[i].split(' ')[0]:
                    shape = int(orders[i].split(' ')[0][-1])
                    title = r'$\bar{u}_v~%s$' % (orders[i].split(' ')[1])
                elif 'u' in orders[i].split(' ')[0]:
                    shape = int(orders[i].split(' ')[0][-1])
                    title = r'$u~%s$' % (orders[i].split(' ')[1])
                elif 'g' in orders[i].split(' ')[0]:
                    shape = int(orders[i].split(' ')[0][-1])
                    title = r'$g~%s$' % (orders[i].split(' ')[1])
                elif 'dv' in orders[i].split(' ')[0]:
                    shape = int(orders[i].split(' ')[0][-1])
                    title = r'$d_v^{\left( %d \right)}~%s$' % (shape, orders[i].split(' ')[1])
                elif orders[i].split(' ')[0][1] == 'b':
                    shape = int(orders[i].split(' ')[0][-1])
                    title = r'$\overline{%s}^{\left( %d \right)}~%s$' % (orders[i].split(' ')[0][0], shape, orders[i].split(' ')[1])
                elif orders[i].split(' ')[0][1] == 'p':
                    shape = int(orders[i].split(' ')[0][-1])
                    title = r'$%s_+^{\left( %d \right)}~%s$' % (orders[i].split(' ')[0][0], shape, orders[i].split(' ')[1])
                elif orders[i].split(' ')[0][1] == 'm':
                    shape = int(orders[i].split(' ')[0][-1])
                    title = r'$%s_-^{\left( %d \right)}~%s$' % (orders[i].split(' ')[0][0], shape, orders[i].split(' ')[1])
                elif 'sea' in orders[i].split(' ')[0]:
                    shape = int(orders[i].split(' ')[0][-1])
                    title = r'$\mathrm{sea}^{\left( %d \right)}~%s$' % (shape, orders[i].split(' ')[1])
                else:
                    #shape = int(orders[i].split(' ')[0][-1])
                    shape=''
                    #title = r'$%s^{\left( %d \right)}~%s$' % (orders[i].split(' ')[0][0], shape, orders[i].split(' ')[1])
                    #title = r'$%s^{\left( %d \right)}~%s$' % (orders[i].split(' ')[0][0], shape, orders[i].split(' ')[1])
                    title='lambda'
                ax.set_title(title, size = 30)
                ## if above chunck causes errors, comment above out and use instead
                ## ax.set_title('%s'%(orders[i]), size=20)
                if histogram:
                    for j in range(len(set(clusters))):
                        color  = colors[j]
                        par = [params[i][k] for k in range(len(params[i])) if clusters[k]==j]
                        ax.hist(par,color=color,alpha=0.6,edgecolor='black')
                        ax.axvline(np.mean(par),ymin=0,ymax=1,ls='--',color=color,alpha=0.8)
                        ax.axvspan(np.mean(par)-np.std(par),np.mean(par)+np.std(par),alpha=0.1,color=color)
                else:
                    for j in range(len(set(clusters))):
                        color  = colors[j]
                        par = [params[i][k] for k in range(len(params[i])) if clusters[k] == j]
                        xs = [_ + 1 for _ in range(len(params[i])) if clusters[_] == j]
                        ax.scatter(xs, par, color = color)
                        ax.axhline(np.mean(par), ls = '--', alpha = 0.5, color = color)
                        ax.axhspan(np.mean(par) - np.std(par), np.mean(par) + np.std(par), alpha = 0.1, color = color)

                i_ax += 1
        i_row += 1

    if histogram: filename='%s/gallery/%s-params-hist.png'%(wdir,dist)
    else:         filename='%s/gallery/%s-params.png'%(wdir,dist)
    checkdir('%s/gallery'%wdir)
    py.tight_layout()
    py.savefig(filename)
    print 'saving figure to %s'%filename




