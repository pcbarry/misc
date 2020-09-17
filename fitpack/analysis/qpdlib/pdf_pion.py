import sys,os
import numpy as np
import copy
from subprocess import Popen, PIPE, STDOUT

#--matplotlib
import matplotlib
matplotlib.use('Agg')
#matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
#matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#matplotlib.rc('text',usetex=True)
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
import pylab as py
import matplotlib.gridspec as gridspec

#--from scipy stack 
from scipy.integrate import quad

#--from tools
from tools.tools     import load,save,checkdir,lprint
from tools.config    import conf,load_config

#--from fitlib
from fitlib.resman import RESMAN

#--from local
from analysis.corelib import core
from analysis.corelib import classifier

FLAV=[]
FLAV.append('g')
FLAV.append('valence')
FLAV.append('sea')

def gen_xf(wdir,Q2=None):
    
    print('\ngenerating pdf-pions from %s'%(wdir))

    load_config('%s/input.py'%wdir)
    istep=core.get_istep()
    replicas=core.get_replicas(wdir)
    core.mod_conf(istep,replicas[0]) #--set conf as specified in istep   

    if 'pdf-pion' not in conf['steps'][istep]['active distributions']:
        print('pdf-pion not in active distribution')
        return 

    resman=RESMAN(nworkers=1,parallel=False,datasets=False)
    parman=resman.parman

    jar=load('%s/data/jar-%d.dat'%(wdir,istep))
    jar_replicas=jar['replicas']
    parman.order=jar['order']

    pdf=conf['pdf-pion']

    #--setup kinematics
    X=10**np.linspace(-3,-1,100)
    X=np.append(X,np.linspace(0.1,0.99,100))
    if Q2==None: Q2=conf['Q20']

    #--compute XF for all replicas        
    XF={}
    cnt=0
    for par in jar_replicas:
        cnt+=1
        lprint('%d/%d'%(cnt,len(replicas)))

        ##--filter
        #flag=False
        #params=replica['params'][istep]
        #order=replica['order'][istep]
        #for i in range(len(order)):
        #    if order[i][0]!=1:continue
        #    if order[i][1]!='pdf':continue
        #    #if order[i][2]=='s1 a':
        #    #   if params[i]<-0.9: flag=True
        #if flag: continue
        core.mod_conf(istep,replicas[cnt-1])

        parman.set_new_params(par,initial=True)

        for flav in FLAV:
            if flav not in XF:  XF[flav]=[]

            if   flav=='valence':
                 func=lambda x: pdf.get_xF(x,Q2,'ub') - pdf.get_xF(x,Q2,'u')
            elif flav=='sea': 
                 func=lambda x: pdf.get_xF(x,Q2,'u') 
            else:
                 func=lambda x: pdf.get_xF(x,Q2,flav) 

            XF[flav].append([func(x) for x in X])
    print     
    checkdir('%s/data'%wdir)
    if Q2==conf['Q20']:
        save({'X':X,'Q2':Q2,'XF':XF},'%s/data/pdf-%d.dat'%(wdir,istep))
    else:
        save({'X':X,'Q2': Q2,'XF':XF},'%s/data/pdf-%d-%d.dat'%(wdir,istep,int(Q2)))

def plot_xf(wdir,kc,Q2=None):

    load_config('%s/input.py'%wdir)
    istep=core.get_istep()
    if Q2==None: Q2=conf['Q20']

    if Q2==conf['Q20']: data=load('%s/data/pdf-%d.dat'%(wdir,istep))
    else: data=load('%s/data/pdf-%d-%d.dat'%(wdir,istep,int(Q2)))

    cluster,colors,nc,cluster_order = classifier.get_clusters(wdir,istep,kc) 
    best_cluster=cluster_order[0]
  

    nrows,ncols=1,1
    fig = py.figure(figsize=(ncols*7,nrows*5))
    ax=py.subplot(nrows,ncols,1)

    X=data['X']
    for i in range(len(data['XF']['g'])):
        if cluster[i]!=best_cluster: continue
        f=data['XF']['g'][i]
        ax.plot(X,np.array(f)/10,'r-',alpha=0.1)

    for i in range(len(data['XF']['valence'])):
        if cluster[i]!=best_cluster: continue
        f=data['XF']['valence'][i]
        ax.plot(X,np.array(f),'g-',alpha=0.1)

    for i in range(len(data['XF']['sea'])):
        if cluster[i]!=best_cluster: continue
        f=data['XF']['sea'][i]
        ax.plot(X,np.array(f),'b-',alpha=0.1)

    ax.set_xlim(1e-3,1)
    ax.semilogx()
    ax.set_ylim(0,0.6)

    ax.tick_params(axis='both', which='major', labelsize=20)
    #ax.legend(loc=2,bbox_to_anchor=(0.5,1))
    ax.set_xlabel(r'$x_{\pi}$',size=25)
    ax.set_ylabel(r'$x_{\pi}f(x_{\pi})$',size=25)
    ax.set_yticks([0.1,0.2,0.3,0.4,0.5,0.6])
    ax.set_xticks([0.001,0.01,0.1,1])
    ax.set_xticklabels([r'$0.001$',r'$0.01$',r'$0.1$',r'$1$'])

    py.tight_layout()
    checkdir('%s/gallery'%wdir)
    py.savefig('%s/gallery/pdf-pion-1-%d.png'%(wdir,istep))

def plot_xf(wdir,kc,Q2=None):

    load_config('%s/input.py'%wdir)
    istep=core.get_istep()
    if Q2==None: Q2=conf['Q20']

    if Q2==conf['Q20']: data=load('%s/data/pdf-%d.dat'%(wdir,istep))
    else: data=load('%s/data/pdf-%d-%d.dat'%(wdir,istep,int(Q2)))

    cluster,colors,nc,cluster_order = classifier.get_clusters(wdir,istep,kc) 
    best_cluster=cluster_order[0]
  
    nrows,ncols=3,1
    fig = py.figure(figsize=(ncols*7,nrows*4))

    X=data['X']
    ax1=py.subplot(nrows,ncols,1)
    for i in range(len(data['XF']['g'])):
        if cluster[i]!=best_cluster: continue
        f=data['XF']['g'][i]
        ax1.plot(X,np.array(f)/10,'r-',alpha=0.1)
    ax1.plot(X,np.mean(data['XF']['g'],axis=0)/10,'k-')

    ax2=py.subplot(nrows,ncols,2)
    for i in range(len(data['XF']['valence'])):
        if cluster[i]!=best_cluster: continue
        f=data['XF']['valence'][i]
        ax2.plot(X,np.array(f),'g-',alpha=0.1)
    ax2.plot(X,np.mean(data['XF']['valence'],axis=0),'k-')

    ax3=py.subplot(nrows,ncols,3)
    for i in range(len(data['XF']['sea'])):
        if cluster[i]!=best_cluster: continue
        f=data['XF']['sea'][i]
        ax3.plot(X,np.array(f),'b-',alpha=0.1)
    ax3.plot(X,np.mean(data['XF']['sea'],axis=0),'k-')

    for ax in [ax1,ax2,ax3]:
        ax.set_xlim(1e-3,1)
        ax.semilogx()
        ax.set_ylim(-0.5,0.5)

        ax.tick_params(axis='both', which='major', labelsize=20)
        #ax.legend(loc=2,bbox_to_anchor=(0.5,1))
        #ax.set_yticks([0.1,0.2,0.3,0.4,0.5,0.6])
        ax.set_xticks([0.001,0.01,0.1,1])
        ax.set_yticks([-0.5,0.0,0.5])
        ax.set_xticklabels([r'$0.001$',r'$0.01$',r'$0.1$',r'$1$'])

    ax1.set_ylabel(r'$x_{\pi}f(x_{\pi})$',size=25)
    ax3.set_xlabel(r'$x_{\pi}$',size=25)

    ax1.text(0.1,0.1,r'\boldmath{$\rm glue/10$}', transform=ax1.transAxes,size=25)
    ax2.text(0.1,0.1,r'\boldmath{$\rm valence$}', transform=ax2.transAxes,size=25)
    ax3.text(0.1,0.1,r'\boldmath{$\rm sea$}', transform=ax3.transAxes,size=25)

    py.tight_layout()
    checkdir('%s/gallery'%wdir)
    py.savefig('%s/gallery/pdf-pion-2-%d.png'%(wdir,istep))

def plot_xf(wdir,kc,Q2=None):

    load_config('%s/input.py'%wdir)
    istep=core.get_istep()
    if Q2==None: Q2=conf['Q20']

    if Q2==conf['Q20']: data=load('%s/data/pdf-%d.dat'%(wdir,istep))
    else: data=load('%s/data/pdf-%d-%d.dat'%(wdir,istep,int(Q2)))

    cluster,colors,nc,cluster_order = classifier.get_clusters(wdir,istep,kc) 
    best_cluster=cluster_order[0]
  

    nrows,ncols=1,1
    fig = py.figure(figsize=(ncols*7,nrows*5))
    ax=py.subplot(nrows,ncols,1)

    X=data['X']
    for i in range(len(data['XF']['g'])):
        if cluster[i]!=best_cluster: continue
        f=data['XF']['g'][i]
        ax.plot(X,np.array(f)/10,'r-',alpha=0.1)

    for i in range(len(data['XF']['valence'])):
        if cluster[i]!=best_cluster: continue
        f=data['XF']['valence'][i]
        ax.plot(X,np.array(f),'g-',alpha=0.1)

    for i in range(len(data['XF']['sea'])):
        if cluster[i]!=best_cluster: continue
        f=data['XF']['sea'][i]
        ax.plot(X,np.array(f),'b-',alpha=0.1)

    ax.set_xlim(1e-3,1)

    ax.tick_params(axis='both', which='major', labelsize=20)
    #ax.legend(loc=2,bbox_to_anchor=(0.5,1))
    ax.set_xlabel(r'$x_{\pi}$',size=25)
    ax.set_ylabel(r'$x_{\pi}f(x_{\pi})$',size=25)
    ax.set_yticks([0.1,0.2,0.3,0.4,0.5,0.6])

    py.tight_layout()
    checkdir('%s/gallery'%wdir)
    py.savefig('%s/gallery/pdf-pion-2-%d.png'%(wdir,istep))

def plot_xf(wdir,kc,Q2=None):

    load_config('%s/input.py'%wdir)
    istep=core.get_istep()

    if Q2==None: Q2=conf['Q20']

    if Q2==conf['Q20']: data=load('%s/data/pdf-%d.dat'%(wdir,istep))
    else: data=load('%s/data/pdf-%d-%d.dat'%(wdir,istep,int(Q2)))

    cluster,colors,nc,cluster_order = classifier.get_clusters(wdir,istep,kc) 
    best_cluster=cluster_order[0]
 
    nrows,ncols=3,1
    fig = py.figure(figsize=(ncols*7.2,nrows*4))
    #gs = fig.add_gridspec(nrows=3, ncols=2, left=0.05, right=0.48, wspace=0.05)
    gs = gridspec.GridSpec(ncols=2, nrows=3, figure=fig
         ,left=0.17, right=0.95,top=0.99,bottom=0.05 ,wspace=0,hspace=0.05)
    ax11 = fig.add_subplot(gs[0,0])
    ax12 = fig.add_subplot(gs[0,1])
    ax21 = fig.add_subplot(gs[1,0])
    ax22 = fig.add_subplot(gs[1,1])
    ax31 = fig.add_subplot(gs[2,0])
    ax32 = fig.add_subplot(gs[2,1])

    def plot_pdf(axl,axr,X,Y,ymin,ymax,factor=1,color='r'):
        for i in range(len(Y)):
            #if cluster[i]!=best_cluster: continue
            axl.plot(X,np.array(Y[i])*factor,color=color,ls='-',alpha=0.1)
            axr.plot(X,np.array(Y[i])*factor,color=color,ls='-',alpha=0.1)
        axl.plot(X,np.mean(Y,axis=0)*factor,'k-',lw=3)
        axr.plot(X,np.mean(Y,axis=0)*factor,'k-',lw=3)
        axl.semilogx()
        axl.set_xlim(1e-3,0.1)
        axr.set_xlim(0.1,1)
        axl.set_ylim(ymin,ymax)
        axr.set_ylim(ymin,ymax)
        axr.set_yticklabels([])
        axr.set_yticks([])
        axr.spines['left'].set_visible(False)
        axl.spines['right'].set_visible(False)
        axl.tick_params(axis='both', which='both',labelsize=20,direction='inout',length=5)
        axr.tick_params(axis='both', which='both',labelsize=20,direction='inout',length=5)
        axl.set_xticks([1e-2,1e-1])
        axr.set_xticks([0.4,0.6,0.8])

    X=data['X']
    ymin,ymax=-0.2,0.5

    Y=data['XF']['g']
    plot_pdf(ax11,ax12,X,Y,ymin,ymax,factor=0.1,color='r')
    Y=data['XF']['valence']
    plot_pdf(ax21,ax22,X,Y,ymin,ymax,factor=1.0,color='g')
    Y=data['XF']['sea']
    plot_pdf(ax31,ax32,X,Y,ymin,ymax,factor=1.0,color='b')


    ax22.plot(X,(1-X),'k:',lw=3,label=r'$(1-x)$')
    ax22.plot(X,(1-X)**2,'k--',lw=3,label=r'$(1-x)^2$')
    ax21.plot(X,(1-X),'k:',lw=3,label=r'$(1-x)$')
    ax21.plot(X,(1-X)**2,'k--',lw=3,label=r'$(1-x)^2$')
    ax21.legend(loc=3,fontsize=20)

    for ax in [ax11,ax12,ax21,ax22]:
        ax.set_xticklabels([])

    ax11.set_ylabel(r'$x_{\pi}f(x_{\pi})$',size=30)
    ax11.yaxis.set_label_coords(-0.2, 0.5)
    ax32.set_xlabel(r'$x_{\pi}$',size=30)
    ax32.xaxis.set_label_coords(0.95, -0.05)

    ax31.set_xticklabels([r'$0.01$',r'$0.1$'])

    props = dict(facecolor='white', alpha=1.0)

    ax11.text(0.1,0.8,r'\boldmath{$\rm glue/10$}'
             ,transform=ax11.transAxes,size=25,bbox=props)
    ax21.text(0.1,0.8,r'\boldmath{$\rm valence$}'
             ,transform=ax21.transAxes,size=25,bbox=props)
    ax31.text(0.1,0.8,r'\boldmath{$\rm sea$}'
             ,transform=ax31.transAxes,size=25,bbox=props)

    py.tight_layout()
    checkdir('%s/gallery'%wdir)
    if Q2!=None: 
        py.savefig('%s/gallery/pdf-pion-%d.png'%(wdir,istep))
    else: 
        py.savefig('%s/gallery/pdf-pion-Q20-%d.png'%(wdir,istep))

def plot_xf(wdir,kc,Q2=None):

    load_config('%s/input.py'%wdir)
    istep=core.get_istep()

    if Q2==None: Q2=conf['Q20']

    if Q2==conf['Q20']: data=load('%s/data/pdf-%d.dat'%(wdir,istep))
    else: data=load('%s/data/pdf-%d-%d.dat'%(wdir,istep,int(Q2)))

    cluster,colors,nc,cluster_order = classifier.get_clusters(wdir,istep,kc) 
    best_cluster=cluster_order[0]
 
    nrows,ncols=1,3
    fig = py.figure(figsize=(ncols*7.2,nrows*4))
    #gs = fig.add_gridspec(nrows=3, ncols=2, left=0.05, right=0.48, wspace=0.05)
    gs = gridspec.GridSpec(ncols=6, nrows=1, figure=fig
         ,left=0.1, right=0.95,top=0.95,bottom=0.25 ,wspace=0.0,hspace=0.05)
    ax11 = fig.add_subplot(gs[0,0])
    ax12 = fig.add_subplot(gs[0,1])
    ax21 = fig.add_subplot(gs[0,2])
    ax22 = fig.add_subplot(gs[0,3])
    ax31 = fig.add_subplot(gs[0,4])
    ax32 = fig.add_subplot(gs[0,5])

    def plot_pdf(axl,axr,X,Y,ymin,ymax,factor=1,color='r'):

        for i in range(len(Y)):
            if cluster[i]!=best_cluster: continue
            #color=colors[cluster[i]]
            axl.plot(X,np.array(Y[i])*factor,color=color,ls='-',alpha=0.2)
            axr.plot(X,np.array(Y[i])*factor,color=color,ls='-',alpha=0.2)

        mean=np.mean(Y,axis=0)*factor
        std=np.std(Y,axis=0)*factor

        #axl.fill_between(X,mean-std,mean+std
        #      ,facecolor=None
        #      #,edgecolor='w'
        #      ,color='Y'
        #      ,alpha=1,zorder=0,hatch='...')
        #axr.fill_between(X,mean-std,mean+std
        #      ,facecolor=None
        #      #,edgecolor='Y'
        #      ,color='Y'
        #      ,alpha=1,zorder=0,hatch='...')

        axl.plot(X,mean+std,'k-',lw=3)
        axl.plot(X,mean-std,'k-',lw=3)
        axr.plot(X,mean+std,'k-',lw=3)
        axr.plot(X,mean-std,'k-',lw=3)

        axl.plot(X,mean,'k--',lw=3)
        axr.plot(X,mean,'k--',lw=3)

        axl.semilogx()
        axl.set_xlim(1e-3,0.1)
        axr.set_xlim(0.1,1)
        axl.set_ylim(ymin,ymax)
        axr.set_ylim(ymin,ymax)
        axr.set_yticklabels([])
        axr.set_yticks([])
        axr.spines['left'].set_visible(False)
        axl.spines['right'].set_visible(False)
        axl.tick_params(axis='both', which='both',labelsize=20,direction='inout',length=5)
        axr.tick_params(axis='both', which='both',labelsize=20,direction='inout',length=5)
        axl.set_xticks([1e-2,1e-1])
        axr.set_xticks([0.4,0.6,0.8])

    X=data['X']
    ymin,ymax=0,0.5

    Y=data['XF']['g']
    plot_pdf(ax11,ax12,X,Y,ymin,ymax,factor=0.1,color='r')
    Y=data['XF']['valence']
    plot_pdf(ax21,ax22,X,Y,ymin,ymax,factor=1.0,color='g')
    Y=data['XF']['sea']
    plot_pdf(ax31,ax32,X,Y,ymin,ymax,factor=1.0,color='b')


    ax22.plot(X,(1-X),'m:',lw=3,label=r'$(1-x)$')
    ax22.plot(X,(1-X)**2,'m--',lw=3,label=r'$(1-x)^2$')
    ax21.plot(X,(1-X),'m:',lw=3,label=r'$(1-x)$')
    ax21.plot(X,(1-X)**2,'m--',lw=3,label=r'$(1-x)^2$')
    ax22.legend(loc=3,fontsize=20)

    for ax in [ax21,ax22,ax31,ax32]:
        ax.set_yticklabels([])

    ax11.set_ylabel(r'$x_{\pi}f(x_{\pi})$',size=30)
    ax11.yaxis.set_label_coords(-0.2, 0.5)
    ax32.set_xlabel(r'$x_{\pi}$',size=30)
    ax32.xaxis.set_label_coords(0.95, -0.05)

    ax11.set_xticklabels([r'$0.01$',r'$0.1$'])
    ax21.set_xticklabels([r'$0.01$',r'$0.1$'])
    ax31.set_xticklabels([r'$0.01$',r'$0.1$'])

    props = dict(facecolor='white', alpha=1.0)

    ax11.text(0.1,0.8,r'\boldmath{$\rm glue/10$}'
             ,transform=ax11.transAxes,size=25,bbox=props)
    ax21.text(0.1,0.8,r'\boldmath{$\rm valence$}'
             ,transform=ax21.transAxes,size=25,bbox=props)
    ax31.text(0.1,0.8,r'\boldmath{$\rm sea$}'
             ,transform=ax31.transAxes,size=25,bbox=props)

    py.tight_layout()
    checkdir('%s/gallery'%wdir)
    if Q2!=conf['Q20']: 
        py.savefig('%s/gallery/pdf-pion-%d.png'%(wdir,istep))
    else: 
        py.savefig('%s/gallery/pdf-pion-Q20-%d.png'%(wdir,istep))
    py.close()

def plot_xf_rat(wdir,basedir,basestep,kc,Q2=None):

    load_config('%s/input.py'%wdir)
    istep=core.get_istep()

    if Q2==None: Q2=conf['Q20']

    if Q2==conf['Q20']: data=load('%s/data/pdf-%d.dat'%(wdir,istep))
    else: data=load('%s/data/pdf-%d-%d.dat'%(wdir,istep,int(Q2)))

    if Q2==conf['Q20']: base_data=load('%s/data/pdf-%d.dat'%(basedir,basestep))
    else: base_data=load('%s/data/pdf-%d-%d.dat'%(basedir,basestep,int(Q2)))

    cluster,colors,nc,cluster_order = classifier.get_clusters(wdir,istep,kc) 
    best_cluster=cluster_order[0]
 
    nrows,ncols=1,3
    fig = py.figure(figsize=(ncols*7.2,nrows*4))
    #gs = fig.add_gridspec(nrows=3, ncols=2, left=0.05, right=0.48, wspace=0.05)
    gs = gridspec.GridSpec(ncols=6, nrows=1, figure=fig
         ,left=0.1, right=0.95,top=0.95,bottom=0.25 ,wspace=0.0,hspace=0.05)
    ax11 = fig.add_subplot(gs[0,0])
    ax12 = fig.add_subplot(gs[0,1])
    ax21 = fig.add_subplot(gs[0,2])
    ax22 = fig.add_subplot(gs[0,3])
    ax31 = fig.add_subplot(gs[0,4])
    ax32 = fig.add_subplot(gs[0,5])

    def plot_pdf(axl,axr,X,Y,ymin,ymax,color='r'):

        for i in range(len(Y)):
            #if cluster[i]!=best_cluster: continue
            axl.plot(X,np.array(Y[i]),color=color,ls='-',alpha=0.1)
            axr.plot(X,np.array(Y[i]),color=color,ls='-',alpha=0.1)

        mean=np.mean(Y,axis=0)
        std=np.std(Y,axis=0)

        #axl.fill_between(X,mean-std,mean+std
        #      ,facecolor=None
        #      #,edgecolor='w'
        #      ,color='Y'
        #      ,alpha=1,zorder=0,hatch='...')
        #axr.fill_between(X,mean-std,mean+std
        #      ,facecolor=None
        #      #,edgecolor='Y'
        #      ,color='Y'
        #      ,alpha=1,zorder=0,hatch='...')

        axl.plot(X,mean+std,'k-',lw=3)
        axl.plot(X,mean-std,'k-',lw=3)
        axr.plot(X,mean+std,'k-',lw=3)
        axr.plot(X,mean-std,'k-',lw=3)

        axl.plot(X,mean,'k--',lw=3)
        axr.plot(X,mean,'k--',lw=3)

        axl.axhline(y=1,linestyle=':')
        axr.axhline(y=1,linestyle=':')

        axl.semilogx()
        axl.set_xlim(1e-3,0.1)
        axr.set_xlim(0.1,1)
        axl.set_ylim(ymin,ymax)
        axr.set_ylim(ymin,ymax)
        axr.set_yticklabels([])
        axr.set_yticks([])
        axr.spines['left'].set_visible(False)
        axl.spines['right'].set_visible(False)
        axl.tick_params(axis='both', which='both',labelsize=20,direction='inout',length=5)
        axr.tick_params(axis='both', which='both',labelsize=20,direction='inout',length=5)
        axl.set_xticks([1e-2,1e-1])
        axr.set_xticks([0.4,0.6,0.8])

    X=data['X']
    ymin,ymax=0.8,1.2

    baseY=np.mean(base_data['XF']['g'],axis=0)
    Y=data['XF']['g']/baseY
    plot_pdf(ax11,ax12,X,Y,ymin,ymax,color='r')
    baseY=np.mean(base_data['XF']['valence'],axis=0)
    Y=data['XF']['valence']/baseY
    plot_pdf(ax21,ax22,X,Y,ymin,ymax,color='g')
    baseY=np.mean(base_data['XF']['sea'],axis=0)
    Y=data['XF']['sea']/baseY
    plot_pdf(ax31,ax32,X,Y,ymin,ymax,color='b')


    for ax in [ax21,ax22,ax31,ax32]:
        ax.set_xticklabels([])

    ax11.set_ylabel(r'$x_{\pi}f(x_{\pi})$',size=30)
    ax11.yaxis.set_label_coords(-0.2, 0.5)
    ax32.set_xlabel(r'$x_{\pi}$',size=30)
    ax32.xaxis.set_label_coords(0.95, -0.05)

    ax11.set_xticklabels([r'$0.01$',r'$0.1$'])
    ax21.set_xticklabels([r'$0.01$',r'$0.1$'])
    ax31.set_xticklabels([r'$0.01$',r'$0.1$'])

    props = dict(facecolor='white', alpha=1.0)

    ax11.text(0.1,0.8,r'\boldmath{$\rm glue$}'
             ,transform=ax11.transAxes,size=25,bbox=props)
    ax21.text(0.1,0.8,r'\boldmath{$\rm valence$}'
             ,transform=ax21.transAxes,size=25,bbox=props)
    ax31.text(0.1,0.8,r'\boldmath{$\rm sea$}'
             ,transform=ax31.transAxes,size=25,bbox=props)

    py.tight_layout()
    checkdir('%s/gallery'%wdir)
    if Q2!=conf['Q20']: 
        py.savefig('%s/gallery/pdfrat-pion-%d.png'%(wdir,istep))
    else: 
        py.savefig('%s/gallery/pdfrat-pion-Q20-%d.png'%(wdir,istep))
    py.close()

def plot_xf_error(DYdir,LNdir,pTdir,kc,Q2=None):

    load_config('%s/input.py'%DYdir)
    DYistep=core.get_istep()

    load_config('%s/input.py'%LNdir)
    LNistep=core.get_istep()

    load_config('%s/input.py'%pTdir)
    pTistep=core.get_istep()

    if Q2==None: Q2=conf['Q20']

    if Q2==conf['Q20']: DYdata=load('%s/data/pdf-%d.dat'%(DYdir,DYistep))
    else: DYdata=load('%s/data/pdf-%d-%d.dat'%(DYdir,DYistep,int(Q2)))

    if Q2==conf['Q20']: LNdata=load('%s/data/pdf-%d.dat'%(LNdir,LNistep))
    else: LNdata=load('%s/data/pdf-%d-%d.dat'%(LNdir,LNistep,int(Q2)))

    if Q2==conf['Q20']: pTdata=load('%s/data/pdf-%d.dat'%(pTdir,pTistep))
    else: pTdata=load('%s/data/pdf-%d-%d.dat'%(pTdir,pTistep,int(Q2)))

    #cluster,colors,nc,cluster_order = classifier.get_clusters(wdir,istep,kc) 
    #best_cluster=cluster_order[0]
 
    nrows,ncols=1,3
    fig = py.figure(figsize=(ncols*7.2,nrows*4))
    #gs = fig.add_gridspec(nrows=3, ncols=2, left=0.05, right=0.48, wspace=0.05)
    gs = gridspec.GridSpec(ncols=6, nrows=1, figure=fig
         ,left=0.1, right=0.95,top=0.95,bottom=0.25 ,wspace=0.0,hspace=0.05)
    ax11 = fig.add_subplot(gs[0,0])
    ax12 = fig.add_subplot(gs[0,1])
    ax21 = fig.add_subplot(gs[0,2])
    ax22 = fig.add_subplot(gs[0,3])
    ax31 = fig.add_subplot(gs[0,4])
    ax32 = fig.add_subplot(gs[0,5])

    #def plot_pdf(axl,axr,X,dat,ymin,ymax,color='r'):
    def plot_pdf(axl,axr,X,dat,ymin,ymax):

        mean={}
        std={}
        for k in ['dy','ln','pT']:
            mean[k]=abs(np.mean(dat[k],axis=0))
            std[k] =np.std(dat[k],axis=0)
            #print mean
            toplot=1+std[k]/mean[k]
            #axl.fill_between(X,toplot,y2=1,color=color,alpha=0.1)
            #axr.fill_between(X,toplot,y2=1,color=color,alpha=0.1)
            if k=='dy':
                axl.fill_between(X,toplot,y2=1,alpha=0.2)
                axr.fill_between(X,toplot,y2=1,alpha=0.2,label=k)
            else:
                axl.fill_between(X,toplot,y2=1,alpha=0.2)
                axr.fill_between(X,toplot,y2=1,alpha=0.2,label='+%s'%k)

            #for i in range(len(mean[k])):
            #    #if cluster[i]!=best_cluster: continue
            #    axl.fill_between(X,np.array(std[k][i]/mean[k][i]),y2=1,color=color,alpha=0.1)
            #    axr.fill_between(X,np.array(std[k][i]/mean[k][i]),y2=1,color=color,alpha=0.1)

        #axl.fill_between(X,mean-std,mean+std
        #      ,facecolor=None
        #      #,edgecolor='w'
        #      ,color='Y'
        #      ,alpha=1,zorder=0,hatch='...')
        #axr.fill_between(X,mean-std,mean+std
        #      ,facecolor=None
        #      #,edgecolor='Y'
        #      ,color='Y'
        #      ,alpha=1,zorder=0,hatch='...')

        axl.axhline(y=1,linestyle=':')
        axr.axhline(y=1,linestyle=':')

        axl.semilogx()
        axl.set_xlim(1e-3,0.1)
        axr.set_xlim(0.1,1)
        axl.set_ylim(ymin,ymax)
        axr.set_ylim(ymin,ymax)
        axr.set_yticklabels([])
        axr.set_yticks([])
        axr.spines['left'].set_visible(False)
        axl.spines['right'].set_visible(False)
        axl.tick_params(axis='both', which='both',labelsize=20,direction='inout',length=5)
        axr.tick_params(axis='both', which='both',labelsize=20,direction='inout',length=5)
        axl.set_xticks([1e-2,1e-1])
        axr.set_xticks([0.4,0.6,0.8])
        #axr.legend(prop={"size":20})

    X=DYdata['X']
    #ymin,ymax=1.0,1.2

    dat={'dy':DYdata['XF']['g'],'ln':LNdata['XF']['g'],'pT':pTdata['XF']['g']}
    ymin,ymax=1.0,2.0
    #plot_pdf(ax11,ax12,X,dat,ymin,ymax,color='r')
    plot_pdf(ax11,ax12,X,dat,ymin,ymax)
    dat={'dy':DYdata['XF']['valence'],'ln':LNdata['XF']['valence'],'pT':pTdata['XF']['valence']}
    ymin,ymax=1.0,1.2
    #plot_pdf(ax21,ax22,X,dat,ymin,ymax,color='g')
    plot_pdf(ax21,ax22,X,dat,ymin,ymax)
    dat={'dy':DYdata['XF']['sea'],'ln':LNdata['XF']['sea'],'pT':pTdata['XF']['sea']}
    ymin,ymax=1.0,2.0
    #plot_pdf(ax31,ax32,X,dat,ymin,ymax,color='b')
    plot_pdf(ax31,ax32,X,dat,ymin,ymax)


    for ax in [ax21,ax22,ax31,ax32]:
        ax.set_yticklabels([])

    ax11.set_ylabel(r'$\sqrt{V[x_{\pi}f(x_{\pi})]}/E[x_\pi f(x_\pi)]$',size=15)
    ax11.yaxis.set_label_coords(-0.2, 0.5)
    ax32.set_xlabel(r'$x_{\pi}$',size=30)
    ax32.xaxis.set_label_coords(0.95, -0.05)

    ax11.set_xticklabels([r'$0.01$',r'$0.1$'])
    ax21.set_xticklabels([r'$0.01$',r'$0.1$'])
    ax31.set_xticklabels([r'$0.01$',r'$0.1$'])

    props = dict(facecolor='white', alpha=1.0)

    ax22.legend(prop={"size":20})

    ax11.text(0.1,0.8,r'\boldmath{$\rm glue$}'
             ,transform=ax11.transAxes,size=25,bbox=props)
    ax21.text(0.1,0.8,r'\boldmath{$\rm valence$}'
             ,transform=ax21.transAxes,size=25,bbox=props)
    ax31.text(0.1,0.8,r'\boldmath{$\rm sea$}'
             ,transform=ax31.transAxes,size=25,bbox=props)

    py.tight_layout()
    checkdir('%s/gallery'%DYdir)
    checkdir('%s/gallery'%LNdir)
    checkdir('%s/gallery'%pTdir)
    if Q2!=conf['Q20']: 
        py.savefig('%s/gallery/pdferror-pion-%d.png'%(DYdir,DYistep))
        py.savefig('%s/gallery/pdferror-pion-%d.png'%(LNdir,LNistep))
        py.savefig('%s/gallery/pdferror-pion-%d.png'%(pTdir,pTistep))
    else: 
        py.savefig('%s/gallery/pdferror-pion-Q20-%d.png'%(DYdir,DYistep))
        py.savefig('%s/gallery/pdferror-pion-Q20-%d.png'%(LNdir,LNistep))
        py.savefig('%s/gallery/pdferror-pion-Q20-%d.png'%(pTdir,pTistep))
    py.close()

def plot_xf_centralrat(DYdir,LNdir,pTdir,kc,Q2=None):

    load_config('%s/input.py'%DYdir)
    DYistep=core.get_istep()

    load_config('%s/input.py'%LNdir)
    LNistep=core.get_istep()

    load_config('%s/input.py'%pTdir)
    pTistep=core.get_istep()

    if Q2==None: Q2=conf['Q20']

    if Q2==conf['Q20']: DYdata=load('%s/data/pdf-%d.dat'%(DYdir,DYistep))
    else: DYdata=load('%s/data/pdf-%d-%d.dat'%(DYdir,DYistep,int(Q2)))

    if Q2==conf['Q20']: LNdata=load('%s/data/pdf-%d.dat'%(LNdir,LNistep))
    else: LNdata=load('%s/data/pdf-%d-%d.dat'%(LNdir,LNistep,int(Q2)))

    if Q2==conf['Q20']: pTdata=load('%s/data/pdf-%d.dat'%(pTdir,pTistep))
    else: pTdata=load('%s/data/pdf-%d-%d.dat'%(pTdir,pTistep,int(Q2)))

    #cluster,colors,nc,cluster_order = classifier.get_clusters(wdir,istep,kc) 
    #best_cluster=cluster_order[0]
 
    nrows,ncols=1,3
    fig = py.figure(figsize=(ncols*7.2,nrows*4))
    #gs = fig.add_gridspec(nrows=3, ncols=2, left=0.05, right=0.48, wspace=0.05)
    gs = gridspec.GridSpec(ncols=6, nrows=1, figure=fig
         ,left=0.1, right=0.95,top=0.95,bottom=0.25 ,wspace=0.0,hspace=0.05)
    ax11 = fig.add_subplot(gs[0,0])
    ax12 = fig.add_subplot(gs[0,1])
    ax21 = fig.add_subplot(gs[0,2])
    ax22 = fig.add_subplot(gs[0,3])
    ax31 = fig.add_subplot(gs[0,4])
    ax32 = fig.add_subplot(gs[0,5])

    def plot_pdf(axl,axr,X,dat,ymin,ymax):

        basemean=np.mean(dat['pT'],axis=0)
        mean={}
        std={}
        #label=[r'$x_F=0.6$','0.7','0.8','0.9']
        for k in ['dy','ln','pT']:
        #for i in range(len(['dy','ln','qT','four'])):
            #k=['dy','ln','qT','four'][i]
            mean[k]=abs(np.mean(dat[k],axis=0))
            std[k] =np.std(dat[k],axis=0)
            #print mean
            toplot=mean[k]/basemean
            if k=='dy':
                axl.plot(X,toplot,label=k)
                axr.plot(X,toplot,label=k)
            else:
                axl.plot(X,toplot,label='+%s'%k)
                axr.plot(X,toplot,label='+%s'%k)
            #axl.plot(X,toplot,label=label[i])
            #axr.plot(X,toplot,label=label[i])
            #axl.fill_between(X,toplot,alpha=0.2)
            #axr.fill_between(X,toplot,alpha=0.2,label=k)

            #for i in range(len(mean[k])):
            #    #if cluster[i]!=best_cluster: continue
            #    axl.fill_between(X,np.array(std[k][i]/mean[k][i]),y2=1,color=color,alpha=0.1)
            #    axr.fill_between(X,np.array(std[k][i]/mean[k][i]),y2=1,color=color,alpha=0.1)

        #axl.fill_between(X,mean-std,mean+std
        #      ,facecolor=None
        #      #,edgecolor='w'
        #      ,color='Y'
        #      ,alpha=1,zorder=0,hatch='...')
        #axr.fill_between(X,mean-std,mean+std
        #      ,facecolor=None
        #      #,edgecolor='Y'
        #      ,color='Y'
        #      ,alpha=1,zorder=0,hatch='...')

        #axl.axhline(y=1,linestyle=':')
        #axr.axhline(y=1,linestyle=':')

        axl.semilogx()
        axl.set_xlim(1e-3,0.1)
        axr.set_xlim(0.1,1)
        axl.set_ylim(ymin,ymax)
        axr.set_ylim(ymin,ymax)
        axr.set_yticklabels([])
        axr.set_yticks([])
        axr.spines['left'].set_visible(False)
        axl.spines['right'].set_visible(False)
        axl.tick_params(axis='both',which='both',labelsize=20,direction='inout',length=5)
        axr.tick_params(axis='both', which='both',labelsize=20,direction='inout',length=5)
        axl.set_xticks([1e-2,1e-1])
        axr.set_xticks([0.4,0.6,0.8])
        axr.legend(prop={"size":20})

    X=DYdata['X']
    ymin,ymax=0.5,1.5

    dat={'dy':DYdata['XF']['g'],'ln':LNdata['XF']['g'],'pT':pTdata['XF']['g']}
    plot_pdf(ax11,ax12,X,dat,ymin,ymax)
    dat={'dy':DYdata['XF']['valence'],'ln':LNdata['XF']['valence'],'pT':pTdata['XF']['valence']}
    plot_pdf(ax21,ax22,X,dat,ymin,ymax)
    dat={'dy':DYdata['XF']['sea'],'ln':LNdata['XF']['sea'],'pT':pTdata['XF']['sea']}
    plot_pdf(ax31,ax32,X,dat,ymin,ymax)


    for ax in [ax21,ax22,ax31,ax32]:
        ax.set_yticklabels([])

    ax11.set_ylabel(r'$f(x_{\pi})/f(x_\pi)^{\rm full}$',size=30)
    ax11.yaxis.set_label_coords(-0.2, 0.5)
    ax32.set_xlabel(r'$x_{\pi}$',size=30)
    ax32.xaxis.set_label_coords(0.95, -0.05)

    ax11.set_xticklabels([r'$0.01$',r'$0.1$'])
    ax21.set_xticklabels([r'$0.01$',r'$0.1$'])
    ax31.set_xticklabels([r'$0.01$',r'$0.1$'])

    #ax32.legend(prop={"size":20})
    props = dict(facecolor='white', alpha=1.0)

    ax11.text(0.1,0.8,r'\boldmath{$\rm glue$}'
             ,transform=ax11.transAxes,size=25,bbox=props)
    ax21.text(0.1,0.8,r'\boldmath{$\rm valence$}'
             ,transform=ax21.transAxes,size=25,bbox=props)
    ax31.text(0.1,0.8,r'\boldmath{$\rm sea$}'
             ,transform=ax31.transAxes,size=25,bbox=props)

    py.tight_layout()
    checkdir('%s/gallery'%DYdir)
    checkdir('%s/gallery'%LNdir)
    checkdir('%s/gallery'%pTdir)
    if Q2!=conf['Q20']: 
        py.savefig('%s/gallery/pdfcentralrat-pion-%d.png'%(DYdir,DYistep))
        py.savefig('%s/gallery/pdfcentralrat-pion-%d.png'%(LNdir,LNistep))
        py.savefig('%s/gallery/pdfcentralrat-pion-%d.png'%(pTdir,pTistep))
    else: 
        py.savefig('%s/gallery/pdfcentralrat-pion-Q20-%d.png'%(DYdir,DYistep))
        py.savefig('%s/gallery/pdfcentralrat-pion-Q20-%d.png'%(LNdir,LNistep))
        py.savefig('%s/gallery/pdfcentralrat-pion-Q20-%d.png'%(pTdir,pTistep))
    py.close()

def plot_xf_central(DYdir,LNdir,pTdir,kc,Q2=None):

    load_config('%s/input.py'%DYdir)
    DYistep=core.get_istep()

    load_config('%s/input.py'%LNdir)
    LNistep=core.get_istep()

    load_config('%s/input.py'%pTdir)
    pTistep=core.get_istep()

    if Q2==None: Q2=conf['Q20']

    if Q2==conf['Q20']: DYdata=load('%s/data/pdf-%d.dat'%(DYdir,DYistep))
    else: DYdata=load('%s/data/pdf-%d-%d.dat'%(DYdir,DYistep,int(Q2)))

    if Q2==conf['Q20']: LNdata=load('%s/data/pdf-%d.dat'%(LNdir,LNistep))
    else: LNdata=load('%s/data/pdf-%d-%d.dat'%(LNdir,LNistep,int(Q2)))

    if Q2==conf['Q20']: qTdata=load('%s/data/pdf-%d.dat'%(pTdir,pTistep))
    else: pTdata=load('%s/data/pdf-%d-%d.dat'%(pTdir,pTistep,int(Q2)))

    #cluster,colors,nc,cluster_order = classifier.get_clusters(wdir,istep,kc) 
    #best_cluster=cluster_order[0]
 
    nrows,ncols=1,3
    fig = py.figure(figsize=(ncols*7.2,nrows*4))
    #gs = fig.add_gridspec(nrows=3, ncols=2, left=0.05, right=0.48, wspace=0.05)
    gs = gridspec.GridSpec(ncols=6, nrows=1, figure=fig
         ,left=0.1, right=0.95,top=0.95,bottom=0.25 ,wspace=0.0,hspace=0.05)
    ax11 = fig.add_subplot(gs[0,0])
    ax12 = fig.add_subplot(gs[0,1])
    ax21 = fig.add_subplot(gs[0,2])
    ax22 = fig.add_subplot(gs[0,3])
    ax31 = fig.add_subplot(gs[0,4])
    ax32 = fig.add_subplot(gs[0,5])

    def plot_pdf(axl,axr,X,dat,ymin,ymax,factor=1):

        mean={}
        std={}
        for k in ['dy','ln','pT']:
            mean[k]=abs(np.mean(dat[k],axis=0))
            std[k] =np.std(dat[k],axis=0)
            #print mean
            if k=='dy': 
                axl.plot(X,mean[k]*factor,label=k)
                axr.plot(X,mean[k]*factor,label=k)
            else:
                axl.plot(X,mean[k]*factor,label='+%s'%k)
                axr.plot(X,mean[k]*factor,label='+%s'%k)
            #axl.fill_between(X,toplot,alpha=0.2)
            #axr.fill_between(X,toplot,alpha=0.2,label=k)

            #for i in range(len(mean[k])):
            #    #if cluster[i]!=best_cluster: continue
            #    axl.fill_between(X,np.array(std[k][i]/mean[k][i]),y2=1,color=color,alpha=0.1)
            #    axr.fill_between(X,np.array(std[k][i]/mean[k][i]),y2=1,color=color,alpha=0.1)

        #axl.fill_between(X,mean-std,mean+std
        #      ,facecolor=None
        #      #,edgecolor='w'
        #      ,color='Y'
        #      ,alpha=1,zorder=0,hatch='...')
        #axr.fill_between(X,mean-std,mean+std
        #      ,facecolor=None
        #      #,edgecolor='Y'
        #      ,color='Y'
        #      ,alpha=1,zorder=0,hatch='...')

        #axl.axhline(y=1,linestyle=':')
        #axr.axhline(y=1,linestyle=':')

        axl.semilogx()
        axl.set_xlim(1e-3,0.1)
        axr.set_xlim(0.1,1)
        axl.set_ylim(ymin,ymax)
        axr.set_ylim(ymin,ymax)
        axr.set_yticklabels([])
        axr.set_yticks([])
        axr.spines['left'].set_visible(False)
        axl.spines['right'].set_visible(False)
        axl.tick_params(axis='both',which='both',labelsize=20,direction='inout',length=5)
        axr.tick_params(axis='both', which='both',labelsize=20,direction='inout',length=5)
        axl.set_xticks([1e-2,1e-1])
        axr.set_xticks([0.4,0.6,0.8])
        axr.legend(prop={"size":20})

    X=DYdata['X']
    ymin,ymax=0.0,0.5

    dat={'dy':DYdata['XF']['g'],'ln':LNdata['XF']['g'],'pT':pTdata['XF']['g']}
    plot_pdf(ax11,ax12,X,dat,ymin,ymax,factor=0.1)
    dat={'dy':DYdata['XF']['valence'],'ln':LNdata['XF']['valence'],'pT':pTdata['XF']['valence']}
    plot_pdf(ax21,ax22,X,dat,ymin,ymax)
    dat={'dy':DYdata['XF']['sea'],'ln':LNdata['XF']['sea'],'pT':pTdata['XF']['sea']}
    plot_pdf(ax31,ax32,X,dat,ymin,ymax)


    for ax in [ax21,ax22,ax31,ax32]:
        ax.set_yticklabels([])

    ax11.set_ylabel(r'$x_{\pi}f(x_{\pi})$',size=30)
    ax11.yaxis.set_label_coords(-0.2, 0.5)
    ax32.set_xlabel(r'$x_{\pi}$',size=30)
    ax32.xaxis.set_label_coords(0.95, -0.05)

    ax11.set_xticklabels([r'$0.01$',r'$0.1$'])
    ax21.set_xticklabels([r'$0.01$',r'$0.1$'])
    ax31.set_xticklabels([r'$0.01$',r'$0.1$'])

    props = dict(facecolor='white', alpha=1.0)

    ax11.text(0.1,0.8,r'\boldmath{$\rm glue/10$}'
             ,transform=ax11.transAxes,size=25,bbox=props)
    ax21.text(0.1,0.8,r'\boldmath{$\rm valence$}'
             ,transform=ax21.transAxes,size=25,bbox=props)
    ax31.text(0.1,0.8,r'\boldmath{$\rm sea$}'
             ,transform=ax31.transAxes,size=25,bbox=props)

    py.tight_layout()
    checkdir('%s/gallery'%DYdir)
    checkdir('%s/gallery'%LNdir)
    checkdir('%s/gallery'%pTdir)
    if Q2!=conf['Q20']: 
        py.savefig('%s/gallery/pdfcentral-pion-%d.png'%(DYdir,DYistep))
        py.savefig('%s/gallery/pdfcentral-pion-%d.png'%(LNdir,LNistep))
        py.savefig('%s/gallery/pdfcentral-pion-%d.png'%(pTdir,pTistep))
    else: 
        py.savefig('%s/gallery/pdfcentral-pion-Q20-%d.png'%(DYdir,DYistep))
        py.savefig('%s/gallery/pdfcentral-pion-Q20-%d.png'%(LNdir,LNistep))
        py.savefig('%s/gallery/pdfcentral-pion-Q20-%d.png'%(pTdir,pTistep))
    py.close()

def plot_xf_many(dirs,kc,Q2=None):

    isteps={}
    for i in range(len(dirs)):
        load_config('%s/input.py'%dirs[i])
        isteps[i]=core.get_istep()

    if Q2==None: Q2=conf['Q20']

    data={}
    if Q2==conf['Q20']: 
        for i in range(len(dirs)):
            data[dirs[i]]=load('%s/data/pdf-%d.dat'%(dirs[i],isteps[i]))
    else: 
        for i in range(len(dirs)):
            data[dirs[i]]=load('%s/data/pdf-%d-%d.dat'%(dirs[i],isteps[i],int(Q2)))

    nrows,ncols=1,3
    fig = py.figure(figsize=(ncols*7.2,nrows*4))
    #gs = fig.add_gridspec(nrows=3, ncols=2, left=0.05, right=0.48, wspace=0.05)
    gs = gridspec.GridSpec(ncols=6, nrows=1, figure=fig
         ,left=0.1, right=0.95,top=0.95,bottom=0.25 ,wspace=0.0,hspace=0.05)
    ax11 = fig.add_subplot(gs[0,0])
    ax12 = fig.add_subplot(gs[0,1])
    ax21 = fig.add_subplot(gs[0,2])
    ax22 = fig.add_subplot(gs[0,3])
    ax31 = fig.add_subplot(gs[0,4])
    ax32 = fig.add_subplot(gs[0,5])

    def plot_pdf(axl,axr,X,Y,ymin,ymax,factor=1):

        colors=['r','g','b','c']
        #colors=['r','g','b','c','m','k']
        scales=[r'$Q$',r'$p_T/2$',r'$p_T$',r'$2p_T$']
        #scales=[r'$x_F=0.6$',r'$0.7$',r'$0.8$',r'$0.9$']
        #scales=[r'$q_{T,{\rm min}}=2.7$',r'$q_{T,{\rm min}}=2.5$',r'$q_{T,{\rm min}}=2.2$',r'$q_{T,{\rm min}}=2.0$',r'$q_{T,{\rm min}}=1.7$',r'$q_{T,{\rm min}}=1.5$']
        #scales=['2.7','2.5','2.2','2.0','1.7','1.5']
        #scales=['2.7','2.2','1.7']
        #for i in range(len(Y)):
        #    k=Y[i]
        #    for j in range(len(k)):
        #        pass
        #        #axl.plot(X,np.array(k[j])*factor,color=colors[i],ls='-',alpha=0.1)
        #        #axr.plot(X,np.array(k[j])*factor,color=colors[i],ls='-',alpha=0.1)
        for i in range(len(Y)):
            #if i==1: continue
            k=Y[i]
            mean=np.mean(k,axis=0)*factor
            std=np.std(k,axis=0)*factor
            #if i==6:
            #if i==3:
            #    axl.plot(X,mean-std,'k-',lw=3)
            #    axl.plot(X,mean+std,'k-',lw=3)
            #    axr.plot(X,mean-std,'k-',lw=3)
            #    axr.plot(X,mean+std,'k-',lw=3)
            #else:
            #axl.plot(X,mean+std,color=colors[i],ls='-',lw=3)
            #axl.plot(X,mean-std,color=colors[i],ls='-',lw=3)
            axl.fill_between(X,mean-std,mean+std,color=colors[i],alpha=0.3)
            #axr.plot(X,mean+std,color=colors[i],ls='-',lw=3)
            #axr.plot(X,mean-std,color=colors[i],ls='-',lw=3)
            axr.fill_between(X,mean-std,mean+std,color=colors[i],alpha=0.3,label=scales[i])



        #axl.fill_between(X,mean-std,mean+std
        #      ,facecolor=None
        #      #,edgecolor='w'
        #      ,color='Y'
        #      ,alpha=1,zorder=0,hatch='...')
        #axr.fill_between(X,mean-std,mean+std
        #      ,facecolor=None
        #      #,edgecolor='Y'
        #      ,color='Y'
        #      ,alpha=1,zorder=0,hatch='...')


        axl.semilogx()
        axl.set_xlim(1e-3,0.1)
        axr.set_xlim(0.1,1)
        axl.set_ylim(ymin,ymax)
        axr.set_ylim(ymin,ymax)
        axr.set_yticklabels([])
        axr.set_yticks([])
        axr.spines['left'].set_visible(False)
        axl.spines['right'].set_visible(False)
        axl.tick_params(axis='both', which='both',labelsize=20,direction='inout',length=5)
        axr.tick_params(axis='both', which='both',labelsize=20,direction='inout',length=5)
        axl.set_xticks([1e-2,1e-1])
        axr.set_xticks([0.4,0.6,0.8])
        #axr.legend(prop={"size":20})

    X=data[dirs[0]]['X']
    ymin,ymax=0.0,0.5

    dat=[]
    for i in range(len(data)):
        dat.append(data[dirs[i]]['XF']['g'])
    plot_pdf(ax11,ax12,X,dat,ymin,ymax,factor=0.1)
    dat=[]
    for i in range(len(data)):
        dat.append(data[dirs[i]]['XF']['valence'])
    plot_pdf(ax21,ax22,X,dat,ymin,ymax)
    dat=[]
    for i in range(len(data)):
        dat.append(data[dirs[i]]['XF']['sea'])
    plot_pdf(ax31,ax32,X,dat,ymin,ymax)

    for ax in [ax21,ax22,ax31,ax32]:
        ax.set_yticklabels([])

    ax11.set_ylabel(r'$x_{\pi}f(x_{\pi})$',size=30)
    ax11.yaxis.set_label_coords(-0.2, 0.5)
    ax32.set_xlabel(r'$x_{\pi}$',size=30)
    ax32.xaxis.set_label_coords(0.95, -0.05)

    ax11.set_xticklabels([r'$0.01$',r'$0.1$'])
    ax21.set_xticklabels([r'$0.01$',r'$0.1$'])
    ax31.set_xticklabels([r'$0.01$',r'$0.1$'])

    ax12.legend(prop={"size":20})
    props = dict(facecolor='white', alpha=1.0)

    ax11.text(0.1,0.8,r'\boldmath{$\rm glue/10$}'
             ,transform=ax11.transAxes,size=25,bbox=props)
    ax21.text(0.1,0.8,r'\boldmath{$\rm valence$}'
             ,transform=ax21.transAxes,size=25,bbox=props)
    ax31.text(0.1,0.8,r'\boldmath{$\rm sea$}'
             ,transform=ax31.transAxes,size=25,bbox=props)

    py.tight_layout()
    for i in range(len(dirs)):
        checkdir('%s/gallery'%dirs[i])
        if Q2!=conf['Q20']:
            py.savefig('%s/gallery/pdfmanyscale-pion-%d.png'%(dirs[i],isteps[i]))
        else: 
            py.savefig('%s/gallery/pdfmanyscale-pion-Q20-%d.png'%(dirs[i],isteps[i]))
    py.close()

def gen_moms(wdir,Q2=None):

    print('\ngenerating pion pdf moments from %s'%(wdir))

    load_config('%s/input.py'%wdir)
    istep=core.get_istep()
    replicas=core.get_replicas(wdir)
    core.mod_conf(istep,replicas[0]) #--set conf as specified in istep   

    if 'pdf-pion' not in conf['steps'][istep]['active distributions']:
        print('pdf-pion not in active distribution')
        return 

    resman=RESMAN(nworkers=1,parallel=False,datasets=False)
    parman=resman.parman

    jar=load('%s/data/jar-%d.dat'%(wdir,istep))
    jar_replicas=jar['replicas']
    parman.order=jar['order']

    pdf=conf['pdf-pion']

    #--setup kinematics
    if Q2==None: Q2=conf['Q20']

    #--compute mom for all replicas        
    mom1={}
    mom2={}
    cnt=0
    for par in jar_replicas:
        cnt+=1
        lprint('%d/%d'%(cnt,len(replicas)))
        core.mod_conf(istep,replicas[cnt-1]) #--set conf as specified in istep   
        parman.set_new_params(par,initial=True)

        #print conf['pdf-pion'].params['g1'][0],conf['pdf-pion'].params['g1'][1]

        for flav in FLAV:
            if flav not in mom1:  mom1[flav]=[]
            if flav not in mom2:  mom2[flav]=[]

            if   flav=='valence':
                 func=lambda x: pdf.get_xF(x,Q2,'ub') - pdf.get_xF(x,Q2,'u')
            elif flav=='sea': 
                 func=lambda x: pdf.get_xF(x,Q2,'u') 
            else:
                 func=lambda x: pdf.get_xF(x,Q2,flav) 

            mom1[flav].append(quad(lambda x: func(x),0,1)[0])
            mom2[flav].append(quad(lambda x: x*func(x),0,1)[0])
    print     
    checkdir('%s/data'%wdir)
    
    if Q2==conf['Q20']: fname='%s/data/pdf-moms-%d.dat'%(wdir,istep)
    else: fname='%s/data/pdf-moms-%d-%d.dat'%(wdir,istep,int(Q2))
    save({'Q2':Q2,'mom1':mom1,'mom2':mom2},fname)

def plot_mom1(wdir,kc,Q2=None):

    load_config('%s/input.py'%wdir)
    istep=core.get_istep()
    if Q2==None: Q2=conf['Q20']

    if Q2==conf['Q20']: data=load('%s/data/pdf-moms-%d.dat'%(wdir,istep))
    else: data=load('%s/data/pdf-moms-%d-%d.dat'%(wdir,istep,int(Q2)))

    cluster,colors,nc,cluster_order = classifier.get_clusters(wdir,istep,kc) 
    best_cluster=cluster_order[0]
  

    nrows,ncols=1,1
    fig = py.figure(figsize=(ncols*7,nrows*5))
    ax=py.subplot(nrows,ncols,1)

    #R=(0,3)
    R=(0,1)
    nbins=50
    mom1=data['mom1']['g']
    mom1=[mom1[i] for i in range(len(mom1)) if cluster[i]==best_cluster]
    ax.hist(mom1,color='r',histtype='step',
            range=R,bins=nbins,label=r'$\rm glue$',density=True)

    mom1=data['mom1']['sea']
    mom1=[6*mom1[i] for i in range(len(mom1)) if cluster[i]==best_cluster]
    ax.hist(mom1,color='b',histtype='step',
            range=R,bins=nbins,label=r'$\rm sea~(\rm tot)$',density=True)

    mom1=data['mom1']['valence']
    mom1=[2*mom1[i] for i in range(len(mom1)) if cluster[i]==best_cluster]
    ax.hist(mom1,color='g',histtype='step',
            range=R,bins=nbins,label=r'$\rm valence~(\rm tot)$',density=True)


    #ax.set_xlim(0,1)
    #ax.set_xlim(0,0.7)
    #ax.set_ylim(0,0.4)

    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.legend(loc=1,fontsize=15)
    ax.set_xlabel(r'$\left<x_{\pi}\right>$',size=25)
    ax.set_ylabel(r'$\rm normalized~yield$',size=25)
    #ax.set_xticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1])
    #ax.set_yticks([0.1,0.2,0.3,0.4])

    py.tight_layout()
    checkdir('%s/gallery'%wdir)
    py.savefig('%s/gallery/pdf-mom-pion-%d.png'%(wdir,istep))

    py.close()

def print_mom1(wdir,kc,Q2=None):
    load_config('%s/input.py'%wdir)
    istep=core.get_istep()
    if Q2==None: Q2=conf['Q20']

    if Q2==conf['Q20']: data=load('%s/data/pdf-moms-%d.dat'%(wdir,istep))
    else: data=load('%s/data/pdf-moms-%d-%d.dat'%(wdir,istep,int(Q2)))

    cluster,colors,nc,cluster_order = classifier.get_clusters(wdir,istep,kc) 
    best_cluster=cluster_order[0]

    best_cluster=cluster_order[0]

    tab={}

    for flav in ['g','sea','valence']:
        
        if    flav=='sea'     : factor=6
        elif  flav=='valence' : factor=2
        else                  : factor=1

        mom1=data['mom1'][flav]
        mom1=[factor*mom1[i] for i in range(len(mom1)) if cluster[i]==best_cluster]
        tab[flav]={}
        tab[flav]['mean']=np.mean(mom1)
        tab[flav]['std']=np.std(mom1)

    print('\n moments from %s\n'%wdir)

    for flav in ['g','sea','valence']:
        msg='%10s = %10.3e +/- %10.3e'
        msg=msg%(flav,tab[flav]['mean'],tab[flav]['std'])
        print(msg)

    tab2={}

    for flav in ['g','sea','valence']:
        if    flav=='sea'     : factor=6
        elif  flav=='valence' : factor=2
        else                  : factor=1

        mom2=data['mom2'][flav]
        mom2=[factor*mom2[i] for i in range(len(mom2)) if cluster[i]==best_cluster]
        tab2[flav]={}
        tab2[flav]['mean']=np.mean(mom2)
        tab2[flav]['std']=np.std(mom2)

    print('\n second moments from %s\n'%wdir)

    for flav in ['g','sea','valence']:
        msg='%10s = %10.3e +/- %10.3e'
        msg=msg%(flav,tab2[flav]['mean'],tab2[flav]['std'])
        print(msg)

def plot_pi2n(wdir,kc):

    print('\n plotting  pi2n lambda from %s'%(wdir))

    load_config('%s/input.py'%wdir)
    istep=core.get_istep()
    replicas=core.get_replicas(wdir)
    core.mod_conf(istep,replicas[0]) #--set conf as specified in istep   

    if 'pdf-pion' not in conf['steps'][istep]['active distributions']:
        print('pdf-pion not in active distribution')
        return 

    Lambda=[]
    for replica in replicas:
        order  = replica['order'][istep]
        params = replica['params'][istep]
        for i in range(len(order)):
            if order[i][1]=='p->pi,n':
                Lambda.append(params[i])

    print np.mean(Lambda)
    print np.std(Lambda)

    nrows,ncols=1,1
    fig = py.figure(figsize=(ncols*7,nrows*5))
    ax=py.subplot(nrows,ncols,1)

    ax.hist(Lambda,color='r',histtype='step',density=True,bins=100)

    #ax.set_xlim(1,2)

    ax.tick_params(axis='both', which='major', labelsize=20)
    #ax.legend(loc=1,fontsize=15)
    ax.set_xlabel(r'$\lambda$',size=25)
    ax.set_ylabel(r'$\rm normalized~yield$',size=25)
    #ax.set_xticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1])
    #ax.set_yticks([0.1,0.2,0.3,0.4])

    py.tight_layout()
    checkdir('%s/gallery'%wdir)
    py.savefig('%s/gallery/pi2n-%d.png'%(wdir,istep))








