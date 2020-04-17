import sys,os
import numpy as np
import copy
from subprocess import Popen, PIPE, STDOUT
import pandas as pd

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

def plot_E615_Q(wdir):

    nrows,ncols=3,4
    fig = py.figure(figsize=(ncols*5,nrows*7))
            
    load_config('%s/input.py'%wdir)
    istep=core.get_istep()

    pred=load('%s/data/predictions-%d.dat'%(wdir,istep))

    pT      =np.array(pred['reactions']['pion_qT'][10000]['pT'])
    total   =np.array(pred['reactions']['pion_qT'][10000]['prediction-rep'])
    asym    =np.array(pred['reactions']['pion_qT'][10001]['prediction-rep'])

    norm = np.ones(total.shape)
    for i in range(len(pred['order'])):
        if pred['order'][i][2] == 1001:
            for j in range(len(pred['params'])/len(pred['order'])):
                norm[j] = norm[j]*pred['params'][i+j*len(pred['order'])]

    print(norm)
    total   = total/norm
    asym    = asym/norm


    data={}
    data['mean']={}
    data['std']={}
    data['mean']['total']=np.mean(total,axis=0)
    data['mean']['asym']=np.mean(asym,axis=0)
    data['std']['total']=np.std(total,axis=0)
    data['std']['asym']=np.std(asym,axis=0)
    data['value']=np.array(pred['reactions']['pion_qT'][1001]['value'])
    data['alpha']=np.array(pred['reactions']['pion_qT'][1001]['alpha'])
    data['pT'] = np.array(pred['reactions']['pion_qT'][1001]['pT'])
    data['Q']=np.array(pred['reactions']['pion_qT'][1001]['Q'])
    Qbins=[]
    #adjust bins
    #Qbins.append([4.05,4.5])
    #Qbins.append([4.5,4.95])
    #Qbins.append([4.95,5.4])
    #Qbins.append([5.4,5.85])
    Qbins.append([5.85,6.75])
    #Qbins.append([6.75,7.65])
    #Qbins=Qbins[::-1]

    f=1 
    for i in range(len(Qbins)):
      ax=py.subplot(nrows,ncols,1)
      X = []
      Y = []
      YM = []
      YP = []
      A = []
      AP = []
      AM = []
      dat = []
      pt = []
      alph = []
      Qavg = (Qbins[i][0]+Qbins[i][1])/2
      for ii in range(len(pT)):
        X.append(pT[ii]/Qavg)
        y = data['mean']['total'][ii]
        Y.append(data['mean']['total'][ii])
        YP.append(y+data['std']['total'][ii])
        YM.append(y-data['std']['total'][ii])
        a = data['mean']['asym'][ii]
        A.append(data['mean']['asym'][ii])
        AP.append(a+data['std']['asym'][ii])
        AM.append(a-data['std']['asym'][ii])
      p,=ax.plot(X,Y,alpha=0.3,label=None)
      ax.fill_between(X,np.abs(YM),np.abs(YP)
        ,color='g',alpha=0.3,label='Thy')
      ax.fill_between(X,np.abs(AM),np.abs(AP)
        ,color='b',alpha=0.3,label='Asy')
      for ii in range(len(data['pT'])):
        if data['Q'][ii]>Qbins[i][0] and data['Q'][ii]<Qbins[i][1]:
          dat.append(data['value'][ii])
          pt.append(data['pT'][ii]/Qavg)
          alph.append(data['alpha'][ii])
      #msg=r'$Q\in[%0.1f,%0.1f]$'%(Qbins[i][0],Qbins[i][1])
          #for k in range(len(data['prediction-rep'])):
          ##for k in range(70):
          #    Yk=d['thy%d'%k].values[ii]
          #    #print k, np.sum(((d['value']-Yk)/d['alpha'])**2)/X.size
          #    ax.plot(X[:-1],(Yk*f**i)[:-1],ls='-'
          #        ,zorder=0
          #        ,color=p.get_color(),alpha=0.3) 

      ax.errorbar(pt,dat,
             alph,zorder=10
              ,marker='o'
              ,ls='none'
              ,markeredgecolor='k'
              ,ecolor='k'
              ,color=p.get_color(),label='data')
      ax.set_title('Q bin: '+str(Qbins[i][0])+'-'+str(Qbins[i][1]))
      #ax.semilogy() 
      #ax.tick_params(axis='both',which='both',direction='in',pad=4,labelsize=18)
      #ax.set_xlabel(r'$pT$',size=20)
      #ax.xaxis.set_label_coords(0.95,-0.03)
      #ax.text(0.06,0.06,r'$d^2\sigma/dQdpT$',transform=ax.transAxes,size=20)
      #ax.text(0.78,0.91,r'$\bf E615 pT$',transform=ax.transAxes,size=17)
      #py.tight_layout()
      ax.legend()
  #        ax.semilogy() 
      ax.tick_params(axis='both',which='both',direction='in',pad=4,labelsize=18)
      ax.set_xlim(None,0.8)
      ax.semilogy()
  #ax.set_ylim(1e-38,1e-20)
  #ax.set_xticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8])
  #ax.set_xticklabels([r'$0$',r'',r'$0.2$',r'',r'$0.4$',r'',r'$0.6$',r'',r''])
  #ax.set_xlabel(r'$pT$',size=20)
      ax.xaxis.set_label_coords(0.95,-0.05)
  #ax.set_ylim(2.5e-1,8e2)
      ax.set_ylabel(r'$d^2\sigma/dQdp_{\rm T}$',size=20)
      ax.set_xlabel(r'$p_{\rm T}/Q$',size=20)
  #        ax.text(0.1,0.23,r'$i=0$',transform=ax.transAxes,size=14)
  #        ax.text(0.5,0.92,r'$i=10$',transform=ax.transAxes,size=14)
      ax.text(0.6,0.9,r'\boldmath{${\rm E615}$}',transform=ax.transAxes,size=17)
    py.tight_layout()
    py.savefig('%s/gallery/dy-pion-E615-QASY-%d.png'%(wdir,istep))

def plot_E615_xF(wdir,istep,data):

    nrows,ncols=1,1
    fig = py.figure(figsize=(ncols*5,nrows*7))
    ax=py.subplot(nrows,ncols,1)            

    Xbins=[]
    #adjust bins
    Xbins.append([0,0.1])
    Xbins.append([0.1,0.2])
    Xbins.append([0.2,0.3])
    Xbins.append([0.3,0.4])
    Xbins.append([0.4,0.5])
    Xbins.append([0.5,0.6])
    Xbins.append([0.7,0.8])
    Xbins.append([0.8,0.9])
    Xbins.append([0.9,1])
    Xbins=Xbins[::-1]

    load_config('%s/input.py'%wdir)
    istep=core.get_istep()
    
    pred=load('%s/data/predictions-%d.dat'%(wdir,istep))

    pT      =np.array(pred['reactions']['pion_qT'][4000]['pT'])
    total   =np.array(pred['reactions']['pion_qT'][4000]['prediction-rep'])
    asym    =np.array(pred['reactions']['pion_qT'][4001]['prediction-rep'])

    data={}
    data['mean']={}
    data['std']={}
    data['mean']['total']=np.mean(total,axis=0)
    data['mean']['asym']=np.mean(asym,axis=0)
    data['std']['total']=np.std(total,axis=0)
    data['std']['asym']=np.std(asym,axis=0)
    data['value']=np.array(pred['reactions']['pion_qT'][4000]['value'])
    data['alpha']=np.array(pred['reactions']['pion_qT'][4000]['alpha'])
    data['xF']=np.array(pred['reactions']['pion_qT'][4000]['alpha'])

    f=1#40 
    for i in range(len(Xbins)):
        d=tab.query('xF>%f and xF<%f'%(Xbins[i][0],Xbins[i][1]))
        if len(d['xF'])!=0:
            ii=np.argsort(d['pT'].values)
            X=d['pT'].values[ii]
            Y=d['thy'].values[ii]
            YP=Y+d['dthy'].values[ii]
            YM=Y-d['dthy'].values[ii]
            if len(d['pT'].values)>1:
                msg=r'$xF\in[%0.1f,%0.1f]$'%(Xbins[i][0],Xbins[i][1])
                p,=ax.plot(X[:-1],(Y*f**i)[:-1],alpha=0.3,label=None)
                ax.fill_between(X[:-1],(YM*f**i)[:-1],(YP*f**i)[:-1]
                  ,color=p.get_color(),alpha=0.3,label=None)

                ax.errorbar(d['pT'][:-1],(d['value']*f**i)[:-1],
                       (d['alpha']*f**i)[:-1],zorder=10
                        ,marker='o'
                        ,ls='none'
                        ,markeredgecolor='k'
                        ,ecolor='k'
                        ,color=p.get_color(),label=msg)
            #ax.set_title('Q bin: '+str(Qbins[i][0])+'-'+str(Qbins[i][1]))
            #ax.semilogy() 
            #ax.tick_params(axis='both',which='both',direction='in',pad=4,labelsize=18)
            #ax.set_xlabel(r'$pT$',size=20)
            #ax.xaxis.set_label_coords(0.95,-0.03)
            #ax.text(0.06,0.06,r'$d^2\sigma/dQdpT$',transform=ax.transAxes,size=20)
            #ax.text(0.78,0.91,r'$\bf E615 pT$',transform=ax.transAxes,size=17)
            #py.tight_layout()
    ax.legend(loc=4)
    ax.semilogy() 
    ax.tick_params(axis='both',which='both',direction='in',pad=4,labelsize=18)
    ax.set_xlim(None,4.5)
    #ax.set_ylim(1e-38,1e-20)
    #ax.set_xticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8])
    #ax.set_xticklabels([r'$0$',r'',r'$0.2$',r'',r'$0.4$',r'',r'$0.6$',r'',r''])
    #ax.set_xlabel(r'$pT$',size=20)
    ax.xaxis.set_label_coords(0.95,-0.03)
    #ax.set_ylim(2.5e-1,8e2)
    ax.set_ylabel(r'$d^2\sigma/dxFdp_{\rm T}~\times %d^i$'%f,size=20)
    ax.set_xlabel(r'$p_{\rm T}$',size=20)
    ax.text(0.1,0.23,r'$i=0$',transform=ax.transAxes,size=14)
    ax.text(0.5,0.92,r'$i=10$',transform=ax.transAxes,size=14)
    ax.text(0.7,0.9,r'\boldmath{${\rm E615}$}',transform=ax.transAxes,size=17)
    py.tight_layout()
    py.savefig('%s/gallery/dy-pion-E615-xF-%d.png'%(wdir,istep))

def plot_obs(wdir,kc): 

    print('\nplotting DY-pion (qT) from %s'%(wdir))

    load_config('%s/input.py'%wdir)
    istep=core.get_istep()
    replicas=core.get_replicas(wdir)
    core.mod_conf(istep,replicas[0]) #--set conf as specified in istep   

    predictions=load('%s/data/predictions-%d.dat'%(wdir,istep))
    #print(predictions['reactions'].keys())
    if 'pion_qT' not in predictions['reactions']: return 

    data=predictions['reactions']['pion_qT']
    plot_E615_Q(wdir,istep,data[1001])
    plot_E615_xF(wdir,istep,data[1002])

    return 

