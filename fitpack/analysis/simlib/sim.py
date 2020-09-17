import os,sys
import pandas as pd
import numpy as np
from subprocess import Popen, PIPE
from scipy.integrate import quad

#--matplotlib
import matplotlib
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
#matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
matplotlib.rc('text',usetex=True)
import pylab  as py
from matplotlib.lines import Line2D

#--from analysis
from analysis.corelib import core

#--from tools
from tools           import config
from tools.tools     import load,save,checkdir,lprint
from tools.config    import conf,load_config
from tools.inputmod  import INPUTMOD
from tools.randomstr import id_generator
from tools.config    import load_config,conf

#--from fitlib
from fitlib.resman import RESMAN
from qcdlib.aux import AUX
from qcdlib.alphaS import ALPHAS

def get_xQ2_grid(option=0,rs=None,nx=10,nQ2=10,fnames=None):

    if option==0:

        Q2=[1.30000E+00,1.50159E+00,1.75516E+00,2.07810E+00\
                ,2.49495E+00,3.04086E+00,3.76715E+00,4.50000E+00\
                ,4.75000E+00,6.23113E+00,8.37423E+00,1.15549E+01\
                ,1.64076E+01,2.40380E+01,3.64361E+01,5.73145E+01\
                ,9.38707E+01,1.60654E+02,2.88438E+02,5.45587E+02\
                ,1.09231E+03,2.32646E+03,5.30043E+03,1.29956E+04\
                ,3.45140E+04,1.00000E+05]

        X=[1.00000E-06,1.28121E-06,1.64152E-06,2.10317E-06\
               ,2.69463E-06,3.45242E-06,4.42329E-06,5.66715E-06\
               ,7.26076E-06,9.30241E-06,1.19180E-05,1.52689E-05\
               ,1.95617E-05,2.50609E-05,3.21053E-05,4.11287E-05\
               ,5.26863E-05,6.74889E-05,8.64459E-05,1.10720E-04\
               ,1.41800E-04,1.81585E-04,2.32503E-04,2.97652E-04\
               ,3.80981E-04,4.87518E-04,6.26039E-04,8.00452E-04\
               ,1.02297E-03,1.30657E-03,1.66759E-03,2.12729E-03\
               ,2.71054E-03,3.44865E-03,4.37927E-03,5.54908E-03\
               ,7.01192E-03,8.83064E-03,1.10763E-02,1.38266E-02\
               ,1.71641E-02,2.11717E-02,2.59364E-02,3.15062E-02\
               ,3.79623E-02,4.53425E-02,5.36750E-02,6.29705E-02\
               ,7.32221E-02,8.44039E-02,9.64793E-02,1.09332E-01\
               ,1.23067E-01,1.37507E-01,1.52639E-01,1.68416E-01\
               ,1.84794E-01,2.01731E-01,2.19016E-01,2.36948E-01\
               ,2.55242E-01,2.73927E-01,2.92954E-01,3.12340E-01\
               ,3.32036E-01,3.52019E-01,3.72282E-01,3.92772E-01\
               ,4.13533E-01,4.34326E-01,4.55495E-01,4.76836E-01\
               ,4.98342E-01,5.20006E-01,5.41818E-01,5.63773E-01\
               ,5.85861E-01,6.08077E-01,6.30459E-01,6.52800E-01\
               ,6.75387E-01,6.98063E-01,7.20830E-01,7.43683E-01\
               ,7.66623E-01,7.89636E-01,8.12791E-01,8.35940E-01\
               ,8.59175E-01,8.82485E-01,9.05866E-01,9.29311E-01\
               ,9.52817E-01,9.76387E-01,1.00000E+00]

    if option==1:

        Q2min=1 #--GeV^2
        Q2max=rs**2
        xmin=Q2min/rs**2
        X=10**np.linspace(np.log10(xmin),0,nx)
        Q2=10**np.linspace(np.log10(Q2min),np.log10(Q2max),nQ2)

    if option==3:
         if fnames !=None:
              for fn in fnames:
                 fdir = os.environ['FITPACK']
                 data= pd.read_excel(fdir + 'database/{}/expdata/{}.xlsx'.format(proc,fdata))
                 X = data['X']
                 Q2= data['Q2']
         else:
              print('No data selected to read')
    return X,Q2

def from_exp_gen_idis_xlsx(wdir,proc,fdata,stat_u_zero=False,syst_u_zero=True):
    """
    rs entry needed only for obs = cross section
    """
    checkdir('%s/sim'%wdir)
    fdir = os.environ['FITPACK']
    #-- the kinem. var.
    print fdata
    data= pd.read_excel(fdir + '/database/{}/expdata/{}.xlsx'.format(proc,fdata))
    lena = len(data[['X']])
    data['col']   =   ['JAM4EIC'] * lena
    data['value'] = [0] * lena
    if stat_u_zero: data['stat_u'] = [1e-10] * lena
    if syst_u_zero: data['syst_u'] = [0.0] * lena

    data = data.loc[:,~data.columns.str.startswith('%cor')]

    df=pd.DataFrame(data)
    df.to_excel('{}/sim/{}.xlsx'.format(wdir,fdata), index=False)

def gen_idis_xlsx(wdir,idx,tar,obs,X,Q2,current=None,lb=None,rs=None,fn=None):
    """
    rs entry needed only for obs = cross section
    """
    checkdir('%s/sim'%wdir)

    #-- the kinem. var.
    data={_:[] for _ in ['col','target','X','Q2','obs','value','stat_u','syst_u']}
    if rs!=None: data['RS']=[]
    if current!=None: data['current']=[]
    if lb!=None: data['lepton beam']=[]

    if fn != None:
      X,Q2=get_xQ2_grid(option=3,fnames=fn)
    else:
      X,Q2=get_xQ2_grid()

    for x in X:
        for q2 in Q2:
            data['col'].append('JAM4EIC')
            data['X'].append(x)
            data['Q2'].append(q2)
            data['target'].append(tar)
            data['obs'].append(obs)
            data['value'].append(0.0)
            data['stat_u'].append(1e-10)
            data['syst_u'].append(0.0)
            if rs!=None: data['RS'].append(rs)
            if current!=None: data['current'].append(current)
            if lb!=None: data['lepton beam'].append(lb)
    df=pd.DataFrame(data)
    df.to_excel('%s/sim/%d.xlsx'%(wdir,idx), index=False)


def gen_conf_idis(wdir,idxs):

    load_config('%s/input.py'%wdir)
    istep=core.get_istep()
    conf['steps'][istep]['datasets']['idis']=[]
    conf['datasets']['idis']['filters']=[]

    for idx in idxs:
        conf['datasets']['idis']['xlsx'][idx]='./%s/sim/%d.xlsx'%(wdir,idx)
        conf['steps'][istep]['datasets']['idis'].append(idx)

    return conf

def update_idis_tabs(wdir,tabs):

    load_config('%s/input.py'%wdir)
    istep=core.get_istep()
    data=load('%s/data/predictions-%d-sim.dat'%(wdir,istep))

    blist=[]
    blist.append('thy')
    blist.append('shift')
    blist.append('residuals')
    blist.append('prediction')
    blist.append('N')
    blist.append('Shift')
    blist.append('W2')
    blist.append('alpha')
    blist.append('residuals-rep')
    blist.append('r-residuals')


    for _ in data['reactions']:
        for idx in data['reactions'][_]:
            print('update %s %d ...'%(_,idx))
            tab=data['reactions'][_][idx]
            for k in blist: del tab[k]

            tab['value']=np.mean(tab['prediction-rep'],axis=0)
            del tab['prediction-rep']
            df=pd.DataFrame(tab)
            df.to_excel('%s/sim/%s.xlsx'%(wdir,idx), index=False)

def get_tables(wdir,idxs):
    table={}
    for idx in idxs:
        tab=pd.read_excel('%s/sim/%d.xlsx'%(wdir,idx))
        tab=tab.to_dict(orient='list')
        table[idx]=tab['value']

    X,Q2=get_xQ2_grid(option=0)
    return X,Q2,table

def gen_lhapdf_dat_file(wdir,idxs,setlabel=0,dirname='JAMSIM'):

    #--get tables
    X,Q2,table=get_tables(wdir,idxs)

    #--start lhapdf file
    lines=[]
    lines.append('PdfType: replica')
    lines.append('Format: lhagrid1')
    lines.append('---')
    line=''
    for _ in X: line+=('%10.5e '%_).upper()
    lines.append(line)
    line=''
    for _ in Q2: line+=('%10.5e '%_**0.5).upper()
    lines.append(line)
    flavs=''
    for idx in idxs: flavs+='%d '%idx
    lines.append(flavs)

    nx=len(X)
    nQ2=len(Q2)
    for i in range(nx*nQ2):
        line=''
        for iflav in idxs:
            line+=('%10.5e '%table[iflav][i]).upper()
        lines.append(line)
    lines.append('---')
    lines=[l+'\n' for l in lines]
    idx=str(setlabel).zfill(4)
    checkdir('%s/lhapdf/%s'%(wdir,dirname))
    tab=open('%s/lhapdf/%s/%s_%s.dat'%(wdir,dirname,dirname,idx),'w')
    tab.writelines(lines)
    tab.close()

def gen_lhapdf_info_file(wdir,idxs,info,dirname='JAMSIM'):

    #--get tables
    X,Q2,table=get_tables(wdir,idxs)

    #--get flav string
    flavs=''
    for idx in idxs: flavs+='%d, '%idx
    flavs=flavs.rstrip(',')

    #--kinematic limits
    xmin=X[0]
    xmax=X[-1]
    Qmin=Q2[0]**0.5
    Qmax=Q2[-1]**0.5

    #--qcd params
    load_config('%s/input.py'%wdir)
    RESMAN(nworkers=1,parallel=False,datasets=False)
    aS=[conf['alphaS'].get_alphaS(_) for _ in Q2]
    mZ=conf['aux'].mZ
    mb=conf['aux'].mb
    mc=conf['aux'].mc
    alphaSMZ=conf['aux'].alphaSMZ

    #--begin lhapdf info file
    lines=[]
    lines.append('SetDesc:         "<description>"')
    lines.append('SetIndex:        <index>')
    lines.append('Authors:         <authors>')
    lines.append('Reference:       <reference>')
    lines.append('Format:          lhagrid1')
    lines.append('DataVersion:     1')
    lines.append('NumMembers:      1')
    lines.append('Particle:        <particle>')
    lines.append('Flavors:         [%s]'%flavs)
    lines.append('OrderQCD:        1')
    lines.append('FlavorScheme:    <flav scheme>')
    lines.append('NumFlavors:      %d'%len(idxs))
    lines.append('ErrorType:       no error')
    lines.append('XMin:            %0.2e'%xmin)
    lines.append('XMax:            %0.2e'%xmax)
    lines.append('QMin:            %0.2e'%Qmin)
    lines.append('QMax:            %0.2e'%Qmax)
    lines.append('MZ:              %f'%mZ)
    lines.append('MUp:             0.0')
    lines.append('MDown:           0.0')
    lines.append('MStrange:        0.0')
    lines.append('MCharm:          %f'%mc)
    lines.append('MBottom:         %f'%mb)
    lines.append('MTop:            180.0')
    lines.append('AlphaS_MZ:       %f'%alphaSMZ)
    lines.append('AlphaS_OrderQCD: 1')
    lines.append('AlphaS_Type:     ipol')
    line='AlphaS_Qs: ['
    for _ in Q2: line+=('%10.5e, '%_**0.5).upper()
    line=line.rstrip(',')+']'
    lines.append(line)
    line='AlphaS_Vals: ['
    for _ in aS: line+=('%10.5e, '%_).upper()
    line=line.rstrip(',')+']'
    lines.append(line)
    lines.append('AlphaS_Lambda4: 0')
    lines.append('AlphaS_Lambda5: 0')

    for i in range(len(lines)):
        for _ in info:
            lines[i]=lines[i].replace(_,info[_])


    lines=[l+'\n' for l in lines]
    checkdir('%s/lhapdf/%s'%(wdir,dirname))
    tab=open('%s/lhapdf/%s/%s.info'%(wdir,dirname,dirname),'w')
    tab.writelines(lines)
    tab.close()

def get_folder_names(wdir):
    process = Popen(['ls','%s/lhapdf'%wdir], stdout=PIPE, stderr=PIPE,universal_newlines=True)
    stdout, stderr = process.communicate()
    names=stdout.split()
    return names

def get_lhapdf_dir():
    try:
        process = Popen(['lhapdf-config','--datadir'], stdout=PIPE, stderr=PIPE,universal_newlines=True)
        stdout, stderr = process.communicate()
        return stdout.strip()
    except:
        sys.exit("lhapdf-config not available")

def benchmark(wdir,idxs,dirname,fmt='pdf'):
    import lhapdf

    X,Q2=get_xQ2_grid(option=0)

    central=lhapdf.mkPDF(dirname,0)

    for idx in idxs:

        nrows,ncols=1,1
        fig = py.figure(figsize=(ncols*4,nrows*3))
        ax=py.subplot(nrows,ncols,1)

        tab=pd.read_excel('%s/sim/%d.xlsx'%(wdir,idx))
        for q2 in Q2:
            bin=tab.query('Q2==%f'%q2)
            X=bin.X.values
            interp = np.array([central.xfxQ2(idx,x,q2) for x in X])
            true   = bin.value.values
            ax.plot(X,np.abs(interp-true)/true*100)

        ax.text(0.1,0.9,r'$%d$'%idx,transform=ax.transAxes)
        ax.set_ylabel(r'$\rm relative~error~(\%)$')
        ax.set_xlabel(r'$x$')
        ax.semilogy()
        ax.semilogx()

        checkdir('%s/gallery'%wdir)
        py.tight_layout()
        py.savefig('%s/gallery/benchmark-%s-%d.%s'%(wdir,dirname,idx,fmt))

#--simualtion of statistical uncertainties

def get_xQ2_bins(rs=None,nx=10,nQ2=10,W2cut=10.,sign=1):

    q2min=1 #--GeV^2
    q2max=rs**2
    xmin=q2min/rs**2
    X=10**np.linspace(np.log10(xmin),0,nx)
    Q2=10**np.linspace(np.log10(q2min),np.log10(q2max),nQ2)

    Xmin,Xmax=X[:-1],X[1:]
    Q2min,Q2max=Q2[:-1],Q2[1:]

    bins=[]

    for ix in range(len(Xmin)):
        for iq2 in range(len(Q2min)):

            #--cut unphysical region
            q2 = Q2max[iq2]
            x  = Xmax[ix]
            q2max=rs**2*x
            if q2>q2max:  continue

            #--impose W2>W2cut
            Q2 = Q2min[iq2]
            x  = Xmax[ix]
            W2 = 0.938**2+ Q2/x*(1-x)
            if W2<W2cut: continue

            b={}
            b['xmin']   = Xmin[ix]
            b['xmax']   = Xmax[ix]
            b['Q2min']  = Q2min[iq2]
            b['Q2max']  = Q2max[iq2]
            b['rs']     = rs
            b['sign']   = sign
            bins.append(b)

    return bins

def plot_kinematics(wdir,data,fmt='pdf'):
    nrows,ncols=1,1
    fig = py.figure(figsize=(ncols*4,nrows*3))
    ax=py.subplot(nrows,ncols,1)

    for b in data:
        X =[b['xmin'],b['xmax'],b['xmax'],b['xmin'],b['xmin']]
        Q2=[b['Q2min'],b['Q2min'],b['Q2max'],b['Q2max'],b['Q2min']]
        ax.plot(X,Q2,'k-')

    ax.set_ylabel(r'$Q^2$')
    ax.set_xlabel(r'$x$')
    ax.semilogy()
    ax.semilogx()

    py.tight_layout()
    checkdir('%s/gallery'%wdir)
    fname='%s/gallery/bins.%s'%(wdir,fmt)
    print fname
    py.savefig(fname)

def get_sigma_dxdy(rs,x,Q2,F2,FL,F3,sign):
    """
    positron sign =-1
    electron sign = 1
    """
    global par
    y=Q2/(rs**2-par.M2)/x
    M2=par.M2
    YP=1+(1-y)**2
    YM=1-(1-y)**2
    thy=2*np.pi*par.alfa**2/x/y/Q2
    thy*=(YP+2*x**2*y**2*M2/Q2)*F2-y**2*FL+sign*YM*x*F3
    return thy

def get_sigma_dxdQ2(rs,x,Q2,F2,FL,F3,sign):
    global par
    dsig=get_sigma_dxdy(rs,x,Q2,F2,FL,F3,sign)
    return dsig/(rs**2-par.M2)/x

def integrad(stf,rs,x,Q2,sign,iF2,iFL,iF3,iw):
    F2=stf.xfxQ2(iF2,x,Q2)
    FL=stf.xfxQ2(iFL,x,Q2)
    F3=0
    if   iw==0: factor=1
    elif iw==1: factor=x
    elif iw==2: factor=Q2
    return get_sigma_dxdQ2(rs,x,Q2,F2,FL,F3,sign)*factor

def integrate(stf,sign,rs,xmin,xmax,q2min,q2max,iF2,iFL,iF3,iw):
    integrand2d=lambda x,q2:integrad(stf,rs,x,q2,sign,iF2,iFL,iF3,iw)
    integrand1d=lambda q2:quad(lambda x:integrand2d(x,q2),xmin,xmax)[0]
    integral=quad(integrand1d,q2min,q2max)[0]
    return integral

def get_sigtot(wdir,dirname,iset,iF2,iFL,iF3,data):
    import lhapdf

    global par
    par=AUX()
    stf=lhapdf.mkPDF(dirname,iset)

    for bin in data:
        xmin =bin['xmin']
        xmax =bin['xmax']
        q2min=bin['Q2min']
        q2max=bin['Q2max']
        sign =bin['sign']
        rs   =bin['rs']
        bin['sigtot']=integrate(stf,sign,rs,xmin,xmax,q2min,q2max,iF2,iFL,iF3,0)
        bin['<x>']=integrate(stf,sign,rs,xmin,xmax,q2min,q2max,iF2,iFL,iF3,1)/bin['sigtot']
        bin['<Q2>']=integrate(stf,sign,rs,xmin,xmax,q2min,q2max,iF2,iFL,iF3,2)/bin['sigtot']

    return data

def red_integrad(stf,rs,x,Q2,sign,iF2,iFL,iF3,iw,current):
    F2=stf.xfxQ2(iF2,x,Q2)
    FL=stf.xfxQ2(iFL,x,Q2)
    F3=0
    aS = ALPHAS()
    M_W = 80.385
    G_F = 1.16638e-5 #value taken from HERAfitter
    y = np.log(np.sqrt(Q2)/(x * rs))
    Yp = 1 + (1 - y)**2
    if   iw==0: factor=1
    elif iw==1: factor=x
    elif iw==2: factor=Q
    if current == 'NC':
        red_factor = (Q2**2 * x)/(2*np.pi*aS**2 * Yp)
    elif current == 'CC':
        red_factor =  (2*np.pi*x/G_F**2)/((M_W**2+Q2)/M_W**2)
    else:
        print('current not specified when calculating sigma reduced')
    return get_sigma_dxdQ2(rs,x,Q2,F2,FL,F3,sign)*factor*red_factor

def red_integrate(stf,sign,rs,xmin,xmax,q2min,q2max,iF2,iFL,iF3,iw,current):
    integrand2d=lambda x,q2:red_integrad(stf,rs,x,q2,sign,iF2,iFL,iF3,iw,current)
    integrand1d=lambda q2:quad(lambda x:integrand2d(x,q2),xmin,xmax)[0]
    integral=quad(integrand1d,q2min,q2max)[0]
    return integral

def get_red_sigtot(wdir,dirname,iset,iF2,iFL,iF3,data,current):
    import lhapdf

    global par
    par=AUX()
    stf=lhapdf.mkPDF(dirname,iset)

    for bin in data:
        xmin =bin['xmin']
        xmax =bin['xmax']
        q2min=bin['Q2min']
        q2max=bin['Q2max']
        sign =bin['sign']
        rs   =bin['rs']
        bin['sigtot']=red_integrate(stf,sign,rs,xmin,xmax,q2min,q2max,iF2,iFL,iF3,0,current)
        bin['<x>']=integrate(stf,sign,rs,xmin,xmax,q2min,q2max,iF2,iFL,iF3,1)/bin['sigtot']
        bin['<Q2>']=integrate(stf,sign,rs,xmin,xmax,q2min,q2max,iF2,iFL,iF3,2)/bin['sigtot']

    return data


def convert_lum(lum):
    one=0.3893793721  #--GeV2 mbarn from PDG
    lum,units=lum.split(':')
    lum=float(lum)
    units=units.strip()
    if units=='fb-1':   return lum*one*1e12
    else:               sys.exit('units not convertible!')

def gen_statistical_uncertainties(wdir,dirname,iset,iF2,iFL,iF3,data,lum):
    """
    ATTENTION: sigtot is in 1/GeV^2
    lum must be in GeV^2
    """
    import lhapdf
    stf=lhapdf.mkPDF(dirname,iset)

    global par
    par=AUX()

    for bin in data:
        N = bin['sigtot']*lum
        ps=(bin['xmax']-bin['xmin'])*(bin['Q2max']-bin['Q2min'])
        bin['<value>']=N/ps/lum
        bin['stat_u']=np.sqrt(N)/ps/lum

        bin['xmid']=0.5*(bin['xmax']  + bin['xmin'])
        bin['Q2mid']=0.5*(bin['Q2max']  + bin['Q2min'])
        x =bin['<x>']
        Q2=bin['<Q2>']
        F2=stf.xfxQ2(iF2,x,Q2)
        FL=stf.xfxQ2(iFL,x,Q2)
        F3=0
        rs=bin['rs']
        sign =bin['sign']
        bin['obs']='dsig/dxdQ2'
        bin['value']=get_sigma_dxdQ2(rs,x,Q2,F2,FL,F3,sign)
        bin['x']=x
        bin['Q2']=Q2
        bin['y']=Q2/x/(rs**2-par.M2)

    return data

def gen_sim_idis_xlsx(wdir,dirname,iset,iF2,iFL,iF3,data,idx,lum,current=None,sig_red=False):

    lum=convert_lum(lum)
    if not sig_red:
        data=get_sigtot(wdir,dirname,iset,iF2,iFL,iF3,data)
    else:
        data=get_red_sigtot(wdir,dirname,iset,iF2,iFL,iF3,data,current)
    data=gen_statistical_uncertainties(wdir,dirname,iset,iF2,iFL,iF3,data,lum)

    cols=[]
    cols.append('col')
    cols.append('X')
    cols.append('Xmid')
    cols.append('Q2')
    cols.append('Q2mid')
    cols.append('RS')
    cols.append('Y')
    cols.append('obs')
    cols.append('value')
    cols.append('stat_u')
    cols.append('<value>')

    table={_:[] for _ in cols}
    for bin in data:
        table['col'].append(dirname)
        table['X'].append(bin['x'])
        table['Xmid'].append(bin['xmid'])
        table['Q2'].append(bin['Q2'])
        table['Q2mid'].append(bin['Q2mid'])
        table['RS'].append(bin['rs'])
        table['Y'].append(bin['y'])
        table['obs'].append(bin['obs'])
        table['value'].append(bin['value'])
        table['stat_u'].append(bin['stat_u'])
        table['<value>'].append(bin['<value>'])

    df=pd.DataFrame(table)
    df.to_excel('%s/sim/%d.xlsx'%(wdir,idx), index=False,columns=cols)
