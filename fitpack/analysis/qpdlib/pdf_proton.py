import sys, os
import numpy as np
import copy
from subprocess import Popen, PIPE, STDOUT

## from scipy stack
from scipy.integrate import quad

## from fitpack tools
from tools.tools     import load, save, checkdir, lprint
from tools.config    import conf, load_config

## from fitpack fitlib
from fitlib.resman import RESMAN

## from fitpack analysis
from analysis.corelib import core
from analysis.corelib import classifier

def gen_xf(wdir, flavors = ['g', 'u', 'ub', 'd', 'db', 's', 'sb'], Q2 = None):
    load_config('%s/input.py' % wdir)
    istep = core.get_istep()

    replicas = core.get_replicas(wdir)
    ## 'conf' will be modified for each replica individually later in the loop over 'replicas'
    ## the reason for doing this is that 'fix parameters' has to be set correctly for each replica

    if 'pdf' not in conf['steps'][istep]['active distributions']:
        print('pdf-proton not in active distribution')
        return

    resman = RESMAN(nworkers = 1, parallel = False, datasets = False)
    parman = resman.parman
    parman.order = replicas[0]['order'][istep]
    ## make sure 'parman' uses the same order for active distributions as all the replicas do
    # print parman.order

    pdf = conf['pdf']

    ## setup kinematics
    X = 10.0 ** np.linspace(-3, -1, 100)
    X = np.append(X, np.linspace(0.1, 0.99, 100))
    if Q2 == None: Q2 = conf['Q20']
    print('\ngenerating pdf-proton from %s at Q2 = %f' % (wdir, Q2))

    ## compute XF for all replicas
    XF = {}
    n_replicas = len(replicas)
    for i in range(n_replicas):
        lprint('%d/%d' % (i + 1, n_replicas))

        ## filter
        #flag=False
        #params=replica['params'][istep]
        #order=replica['order'][istep]
        #for i in range(len(order)):
        #    if order[i][0]!=1:continue
        #    if order[i][1]!='pdf':continue
        #    #if order[i][2]=='s1 a':
        #    #   if params[i]<-0.9: flag=True
        #if flag: continue

        core.mod_conf(istep, replicas[i])
        parman.set_new_params(replicas[i]['params'][istep], initial = True)

        #print conf['pdf-pion'].params['g1'][0],conf['pdf-pion'].params['g1'][1]

        for flavor in flavors:
            if flavor not in XF: XF[flavor] = []
            if flavor == 'rs':
                func = lambda x: (pdf.get_xF(x, Q2, 's') + pdf.get_xF(x, Q2, 'sb')) / (pdf.get_xF(x, Q2, 'db') + pdf.get_xF(x, Q2, 'ub'))
            elif flavor == 'uv':
                func = lambda x: pdf.get_xF(x, Q2, 'u') - pdf.get_xF(x, Q2, 'ub')
            elif flavor == 'dv':
                func = lambda x: pdf.get_xF(x, Q2, 'd') - pdf.get_xF(x, Q2, 'db')
            elif flavor == 'd/u':
                func = lambda x: pdf.get_xF(x, Q2, 'd') / pdf.get_xF(x, Q2, 'u')
            elif flavor == 'db+ub':
                func = lambda x: pdf.get_xF(x, Q2, 'db') + pdf.get_xF(x, Q2, 'ub')
            elif flavor == 'db-ub':
                func = lambda x: pdf.get_xF(x, Q2, 'db') - pdf.get_xF(x, Q2, 'ub')
            elif flavor == 's+sb':
                func = lambda x: pdf.get_xF(x, Q2, 's') + pdf.get_xF(x, Q2, 'sb')
            else:
                func = lambda x: pdf.get_xF(x, Q2, flavor)

            XF[flavor].append([func(x) for x in X])
    print
    checkdir('%s/data' % wdir)
    if Q2 == conf['Q20']:
        save({'X': X, 'Q2': Q2, 'XF': XF}, '%s/data/pdf-%d.dat' % (wdir, istep))
    else:
        save({'X': X, 'Q2': Q2, 'XF': XF}, '%s/data/pdf-%d-%f.dat' % (wdir, istep, Q2))

def get_ns(wdir, Q2 = None):
    load_config('%s/input.py' % wdir)
    istep = core.get_istep()

    replicas = core.get_replicas(wdir)
    ## 'conf' will be modified for each replica individually later in the loop over 'replicas'
    ## the reason for doing this is that 'fix parameters' has to be set correctly for each replica

    if 'pdf' not in conf['steps'][istep]['active distributions']:
        print('pdf-proton not in active distribution')
        return

    resman = RESMAN(nworkers = 1, parallel = False, datasets = False)
    parman = resman.parman
    parman.order = replicas[0]['order'][istep] ## make sure 'parman' uses the same order as all the replicas do
    # print parman.order

    pdf = conf['pdf']

    ## setup kinematics
    if Q2 == None: Q2 = conf['Q20']
    ## get normalization for all flavors of PDF
    print('\ngetting pdf-normalization from %s at Q2 = %f' % (wdir, Q2))

    ns = {'first': {}, 'second': {}}
    n_replicas = len(replicas)
    for i in range(n_replicas):
        lprint('%d/%d' % (i + 1, n_replicas))

        core.mod_conf(istep, replicas[i])
        parman.set_new_params(replicas[i]['params'][istep], initial = True)

        for flavor, value in pdf.params.iteritems():
            if ('1' in flavor) and (flavor not in ns['first']):
                ns['first'][flavor] = []
            elif ('2' in flavor) and (flavor not in ns['second']):
                ns['second'][flavor] = []
            if ('1' in flavor):
                ns['first'][flavor].append(pdf.params[flavor][0])
            elif ('2' in flavor):
                ns['second'][flavor].append(pdf.params[flavor][0])

    print
    checkdir('%s/data' % wdir)
    if Q2 == conf['Q20']:
        save(ns, '%s/data/pdf-normalization-%d.dat' % (wdir, istep))
    else:
        save(ns, '%s/data/pdf-normalization-%d-%f.dat' % (wdir, istep, Q2))

def print_replica_parameters(wdir, istep = 0, i_replica = 1):

    print('\nprinting PDF parameters from %s' % wdir)

    load_config('%s/input.py' % wdir)
    if istep == 0:
        istep = core.get_istep()

    replicas = core.get_replicas(wdir)
    for i in range(len(replicas[i_replica - 1]['order'][istep])):
        print replicas[i_replica - 1]['order'][istep][i], replicas[i_replica - 1]['params'][istep][i]

def print_pdf_parameters(wdir, istep = 0, i_replica = 1):
    ## print parameters in PDF class
    print('\nprinting PDF parameters from %s' % wdir)

    load_config('%s/input.py' % wdir)
    if istep == 0:
        istep = core.get_istep()

    replicas = core.get_replicas(wdir)
    core.mod_conf(istep, replicas[0]) ## set conf as specified in istep

    if 'pdf' not in conf['steps'][istep]['active distributions']:
        print('pdf-proton not in active distribution')
        return

    resman = RESMAN(nworkers = 1, parallel = False, datasets = False)
    parman = resman.parman
    parman.order = replicas[0]['order'][istep] ## make sure 'parman' uses the same order as all the replicas do
    # print parman.order

    pdf = conf['pdf']

    core.mod_conf(istep, replicas[i_replica - 1])
    parman.set_new_params(replicas[i_replica - 1]['params'][istep], initial = True)
    for name, value in pdf.params.iteritems():
        print '%7s: %.5e, %.5e, %.5e, %.5e, %.5e' % (name, value[0], value[1], value[2], value[3], value[4])

if __name__ == '__main__':
    pass
