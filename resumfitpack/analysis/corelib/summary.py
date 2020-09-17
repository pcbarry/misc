import sys,os
import numpy as np
import copy
from subprocess import Popen, PIPE, STDOUT
import pandas as pd
import math

#--matplotlib
import matplotlib
matplotlib.use('Agg')
#matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
#matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#matplotlib.rc('text',usetex=True)
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
import pylab as py


#--from tools
from tools.tools     import load,save,checkdir,lprint
from tools.config    import conf,load_config

#--from fitlib
from fitlib.resman import RESMAN

#--from local
from analysis.corelib import core
from analysis.corelib import classifier

def get_norm(wdir,kc):
    istep=core.get_istep()
    replicas=core.get_replicas(wdir)
    core.mod_conf(istep,replicas[0]) #--set conf as specified in istep

    tab={}
    for replica in replicas:
        order=replica['order'][istep]
        #for _ in order: print(_,istep)
        #sys.exit()
        params=replica['params'][istep]
        for i in range(len(order)):
            if order[i][0]==2:
                reaction=order[i][1]
                idx=order[i][2]
                if reaction not in tab: tab[reaction]={}
                if idx not in tab[reaction]: tab[reaction][idx]=[]
                tab[reaction][idx].append(params[i])

    for reaction in tab:
        for idx in tab[reaction]:
            norm=tab[reaction][idx][:]
            tab[reaction][idx]={}
            tab[reaction][idx]['mean']=np.mean(norm)
            tab[reaction][idx]['std']=np.std(norm)

    return tab

def get_norm_clusters(wdir,kc):
    istep=core.get_istep()
    replicas=core.get_replicas(wdir)
    core.mod_conf(istep,replicas[0]) #--set conf as specified in istep

    cluster,colors,nc,c_order=classifier.get_clusters(wdir,istep,kc)

    tab={}
    for c in c_order:
        tab[c]={}
        for i in range(len(replicas)):
            flag=True
            if cluster[i]==c: flag=False
            if flag: continue
            replica=replicas[i]
            order=replica['order'][istep]
            #for _ in order: print(_,istep)
            #sys.exit()
            params=replica['params'][istep]
            for j in range(len(order)):
                if order[j][0]==2:
                    reaction=order[j][1]
                    idx=order[j][2]
                    if reaction not in tab[c]: tab[c][reaction]={}
                    if idx not in tab[c][reaction]: tab[c][reaction][idx]=[]
                    tab[c][reaction][idx].append(params[j])

        for reaction in tab[c]:
            for idx in tab[c][reaction]:
                norm=tab[c][reaction][idx][:]
                tab[c][reaction][idx]={}
                tab[c][reaction][idx]['mean']=np.mean(norm)
                tab[c][reaction][idx]['std']=np.std(norm)

    return tab

def get_chi2(wdir,kc):

    istep=core.get_istep()
    #replicas=core.get_replicas(wdir)
    #core.mod_conf(istep,replicas[0]) #--set conf as specified in istep

    predictions=load('%s/data/predictions-%d.dat'%(wdir,istep))
    data=predictions['reactions']

    tab={}
    for reaction in data:
        for idx in data[reaction]:
            if reaction not in tab: tab[reaction]={}
            if idx not in tab[reaction]: tab[reaction][idx]={}

            value=data[reaction][idx]['value']
            alpha=data[reaction][idx]['alpha']
            thy=np.mean(data[reaction][idx]['prediction-rep'],axis=0)
            col=data[reaction][idx]['col'][0]
            res=(value-thy)/alpha
            chi2=np.sum(res**2)
            npts=res.size
            chi2_npts=chi2/npts

            tab[reaction][idx]['col']=col
            tab[reaction][idx]['chi2']=chi2
            tab[reaction][idx]['npts']=npts
            tab[reaction][idx]['chi2_npts']=chi2_npts
    return tab

def get_chi2_clusters(wdir,kc):

    istep=core.get_istep()
    #replicas=core.get_replicas(wdir)
    #core.mod_conf(istep,replicas[0]) #--set conf as specified in istep

    predictions=load('%s/data/predictions-%d.dat'%(wdir,istep))
    data=predictions['reactions']

    cluster,colors,nc,c_order=classifier.get_clusters(wdir,istep,kc)

    tab={}
    for c in c_order:
        tab[c]={}
        for reaction in data:
            for idx in data[reaction]:
                if reaction not in tab[c]: tab[c][reaction]={}
                if idx not in tab[c][reaction]: tab[c][reaction][idx]={}

                value=data[reaction][idx]['value']
                alpha=data[reaction][idx]['alpha']
                predictions=[]
                for i in range(len(data[reaction][idx]['prediction-rep'])):
                    if cluster[i]==c: predictions.append(data[reaction][idx]['prediction-rep'][i])
                #thy=np.mean(data[reaction][idx]['prediction-rep'],axis=0)
                thy=np.mean(predictions,axis=0)
                col=data[reaction][idx]['col'][0]
                res=(value-thy)/alpha
                chi2=np.sum(res**2)
                npts=res.size
                chi2_npts=chi2/npts

                tab[c][reaction][idx]['col']=col
                tab[c][reaction][idx]['chi2']=chi2
                tab[c][reaction][idx]['npts']=npts
                tab[c][reaction][idx]['chi2_npts']=chi2_npts
    return tab

def print_summary(wdir,kc):
    load_config('%s/input.py'%wdir)
    istep=core.get_istep()
    summary_file = open("%s/summary.txt" % wdir, "w")
    print('\nsummary of  %s\n'%(wdir))
    load_config('%s/input.py'%wdir)
    norm_tab=get_norm(wdir,kc)
    chi2_tab=get_chi2(wdir,kc)
    msg1 ='%10s '
    msg1+='%10s '
    msg1+='%10s '
    msg1+='%20s '
    msg1+='%5s '
    msg1+='%10s '
    msg1+='%10s '
    msg1+='%10s '
    msg1=msg1%('prediction','reaction','idx','col','npts','chi2','chi2/npts','norm')
    print(msg1)
    print >> summary_file, msg1

    chi2_tot=0
    npts_tot=0
    for reaction in chi2_tab:
        for idx in chi2_tab[reaction]:
            col=chi2_tab[reaction][idx]['col']
            npts=chi2_tab[reaction][idx]['npts']
            chi2=chi2_tab[reaction][idx]['chi2']
            chi2_npts=chi2_tab[reaction][idx]['chi2_npts']
            if reaction in norm_tab:
                if idx in norm_tab[reaction]:
                    norm=norm_tab[reaction][idx]['mean']
                else:
                    norm = None
            else:
                if 'norm' in conf['datasets'][reaction]:
                    if idx in conf['datasets'][reaction]['norm']:
                        norm=conf['datasets'][reaction]['norm'][idx]['value']
                    else:
                        norm = None
                else:
                    norm = None

            prediction=False

            if 'prediction' in conf:
                if reaction in conf['prediction']:
                    if any([idx==_ for _ in conf['prediction'][reaction]]):
                        prediction=True

            if prediction==False:
                chi2_tot+=chi2
                npts_tot+=npts

            msg2 ='%10s '
            msg2+='%10s '
            msg2+='%10d '
            msg2+='%20s '
            msg2+='%5d '
            msg2+='%10.2f '
            msg2+='%10.2f '
            if norm != None:
                msg2+='%10.2f '
                print(msg2 % (prediction,reaction,idx,col,npts,chi2,chi2_npts,norm))
                print >> summary_file, msg2 % (prediction,reaction,idx,col,npts,chi2,chi2_npts,norm)
            else:
                print(msg2 % (prediction,reaction,idx,col,npts,chi2,chi2_npts))
                print >> summary_file, msg2 % (prediction,reaction,idx,col,npts,chi2,chi2_npts)

    chi2_npts_tot=chi2_tot/npts_tot

    print("-"*len(msg1))
    msg3 ='%10s '
    msg3+='%10s '
    msg3+='%10s '
    msg3+='%20s '
    msg3+='%5d '
    msg3+='%10.2f '
    msg3+='%10.2f '
    msg3+='%10s'
    print(msg3 % (' ',' ',' ',' ',npts_tot,chi2_tot,chi2_npts_tot,' '))
    print >> summary_file, msg3 % (' ',' ',' ',' ',npts_tot,chi2_tot,chi2_npts_tot,' ')
    summary_file.close()

def print_summary_clusters(wdir,kc):
    load_config('%s/input.py'%wdir)
    istep=core.get_istep()
    load_config('%s/input.py'%wdir)
    norm_tab=get_norm_clusters(wdir,kc)
    chi2_tab=get_chi2_clusters(wdir,kc)
    cluster,colors,nc,order=classifier.get_clusters(wdir,istep,kc)

    for c in chi2_tab:
        color=colors[c]
        flag=False
        for reaction in chi2_tab[c]:
            for idx in chi2_tab[c][reaction]:
                if math.isnan(chi2_tab[c][reaction][idx]['chi2']): flag=True
        if flag: continue
        summary_file = open("%s/summary-%s.txt" %(wdir,color), "w")
        print('\nsummary of  %s color:%s cluster\n'%(wdir,color))
        msg1 ='%10s '
        msg1+='%10s '
        msg1+='%10s '
        msg1+='%20s '
        msg1+='%5s '
        msg1+='%10s '
        msg1+='%10s '
        msg1+='%10s '
        msg1=msg1%('prediction','reaction','idx','col','npts','chi2','chi2/npts','norm')
        print(msg1)
        print >> summary_file, msg1
        chi2_tot=0
        npts_tot=0
        nexpt   =0
        echi2   =0
        for reaction in chi2_tab[c]:
            for idx in chi2_tab[c][reaction]:
                nexpt+=1
                col=chi2_tab[c][reaction][idx]['col']
                npts=chi2_tab[c][reaction][idx]['npts']
                chi2=chi2_tab[c][reaction][idx]['chi2']
                chi2_npts=chi2_tab[c][reaction][idx]['chi2_npts']
                if reaction in norm_tab[c]:
                    if idx in norm_tab[c][reaction]:
                        norm=norm_tab[c][reaction][idx]['mean']
                    else:
                        norm = None
                else:
                    if 'norm' in conf['datasets'][reaction]:
                        if idx in conf['datasets'][reaction]['norm']:
                            norm=conf['datasets'][reaction]['norm'][idx]['value']
                        else:
                            norm = None
                    else:
                        norm = None

                prediction=False

                if 'prediction' in conf:
                    if reaction in conf['prediction']:
                        if any([idx==_ for _ in conf['prediction'][reaction]]):
                            prediction=True

                if prediction==False:
                    chi2_tot+=chi2
                    npts_tot+=npts
                    echi2   +=chi2/npts

                msg2 ='%10s '
                msg2+='%10s '
                msg2+='%10d '
                msg2+='%20s '
                msg2+='%5d '
                msg2+='%10.2f '
                msg2+='%10.2f '
                if norm != None:
                    msg2+='%10.2f '
                    print(msg2 % (prediction,reaction,idx,col,npts,chi2,chi2_npts,norm))
                    print >> summary_file, msg2 % (prediction,reaction,idx,col,npts,chi2,chi2_npts,norm)
                else:
                    print(msg2 % (prediction,reaction,idx,col,npts,chi2,chi2_npts))
                    print >> summary_file, msg2 % (prediction,reaction,idx,col,npts,chi2,chi2_npts)

        chi2_npts_tot=chi2_tot/npts_tot
        echi2_nexpt =echi2/nexpt

        print("-"*len(msg1))
        msg3 ='%10s '
        msg3+='%10s '
        msg3+='%10s '
        msg3+='%20s '
        msg3+='%5d '
        msg3+='%10.2f '
        msg3+='%10.2f '
        msg3+='%10s'
        print(msg3 % (' ',' ',' ',' ',npts_tot,chi2_tot,chi2_npts_tot,' '))
        print >> summary_file, msg3 % (' ',' ',' ',' ',npts_tot,chi2_tot,chi2_npts_tot,' ')

        print("-"*len(msg1))
        msg4 ='%10s '
        msg4+='%10s '
        msg4+='%10s '
        msg4+='%20s '
        msg4+='%5s '
        msg4+='%10s '
        msg4+='%15s '
        msg4+='%10s'
        msg4=msg4%(' ',' ',' ',' ','nexpt','echi2','echi2/nexpt',' ')
        print(msg4)
        print >> summary_file, msg4 
        msg5 ='%10s '
        msg5+='%10s '
        msg5+='%10s '
        msg5+='%20s '
        msg5+='%5d '
        msg5+='%10.2f '
        msg5+='%10.2f '
        msg5+='%10s'
        print(msg5 % (' ',' ',' ',' ',nexpt,echi2,echi2_nexpt,' '))
        print >> summary_file, msg5 % (' ',' ',' ',' ',nexpt,echi2,echi2_nexpt,' ')
        summary_file.close()





