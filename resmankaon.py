import sys, os
import time
import numpy as np

#--from nuclib
import nuclib

#--from qcdlib
from qcdlib import pdf0, pdf1, pdf2, pdf3, pdf4, pdf10
from qcdlib import ppdf0,ppdf3
from qcdlib import ff0,ff0h,ffsum,ffsum2,ffpion1,ffkaon1
from qcdlib import pdfpion0,pdfpion1,pdfpionsea,pdfkaon
from qcdlib import aux, eweak,ht0,ht1,ht2, alphaS, mellin, offshell0, offshell1, offshell2, offshell3, offshell4

#--from fitlib
from fitlib.parman import PARMAN

#--from tools
from tools.tools    import checkdir
from tools.config   import conf,load_config,options
from tools.parallel import PARALLEL

#--from nuclib
try:
    from nuclib import deuterium
except:
    pass


#--from obslib
try:
    import obslib.idis.residuals
    import obslib.idis.reader
    import obslib.idis.theory
except:
    pass

try:
    import obslib.pidis.residuals
    import obslib.pidis.reader
    import obslib.pidis.theory
except:
    pass

try:
    import obslib.sidis.residuals
    import obslib.sidis.reader
    import obslib.sidis.theory
except:
    pass

try:
    import obslib.psidis.residuals
    import obslib.psidis.reader
    import obslib.psidis.theory
except:
    pass

try:
    import obslib.dy.residuals
    import obslib.dy.reader
    import obslib.dy.dymellin
    import obslib.dy.theory
    import obslib.dy.piontheory
except:
    pass

try:
    import obslib.wzrv.residuals
    import obslib.wzrv.reader
    import obslib.wzrv.wzmellin
    import obslib.wzrv.theory
except:
    pass

try:
    import obslib.wasym.residuals
    import obslib.wasym.reader
    import obslib.wasym.wmellin
    import obslib.wasym.theory
except:
    pass

try:
    import obslib.zrap.residuals
    import obslib.zrap.reader
    import obslib.zrap.zmellin
    import obslib.zrap.theory
except:
    pass

try:
    import obslib.sia.residuals
    import obslib.sia.reader
    import obslib.sia.theory
except:
    pass

try:
    import obslib.jets.residuals
    import obslib.jets.reader
    import obslib.jets.jet_mellin
except:
    pass

try:
    import obslib.pjets.residuals
    import obslib.pjets.reader
    import obslib.pjets.pjet_mellin
except:
    pass

try:
    import obslib.qpdf.residuals
    import obslib.qpdf.reader
    import obslib.qpdf.theory
except:
    pass

try:
    import obslib.pion_qT.residuals
    import obslib.pion_qT.reader
    import obslib.pion_qT.theory
    import obslib.pion_qT.theoryasy
except:
    pass

try:
    import obslib.ln.residuals
    import obslib.ln.reader
    import obslib.ln.theory
except:
    pass

try:
    import obslib.lh.residuals
    import obslib.lh.reader
    import obslib.lh.theory
except:
    pass

class RESMAN:

    def __init__(self,nworkers=2,parallel=True,datasets=True):

        self.setup_core()
        if datasets:
            if 'idis'    in conf['datasets']: self.setup_idis()
            if 'pidis'   in conf['datasets']: self.setup_pidis()
            if 'sidis'   in conf['datasets']: self.setup_sidis()
            if 'psidis'  in conf['datasets']: self.setup_psidis()
            if 'dy'      in conf['datasets']: self.setup_dy()
            if 'wzrv'    in conf['datasets']: self.setup_wzrv()
            if 'wasym'   in conf['datasets']: self.setup_wasym()
            if 'zrap'    in conf['datasets']: self.setup_zrap()
            if 'sia'     in conf['datasets']: self.setup_sia()
            if 'jet'     in conf['datasets']: self.setup_jet()
            if 'pjet'    in conf['datasets']: self.setup_pjet()
            if 'qpdf'    in conf['datasets']: self.setup_qpdf()
            if 'ln'      in conf['datasets']: self.setup_ln()
            if 'dy-pion' in conf['datasets']: self.setup_dy_pion()
            if 'pion_qT' in conf['datasets']: self.setup_pion_qT()
            if 'lh'      in conf['datasets']: self.setup_lh()
        self.parman=PARMAN()

        self.parallel_mode=parallel
        if parallel:
            self.setup_parallel(nworkers)

    def setup_core(self):

        conf['aux'] = aux.AUX()
        if 'mellin npts' in conf: npts = conf['mellin npts']
        else: npts = 4
        conf['mellin']= mellin.MELLIN(npts=npts)
        conf['mellin-pion']=mellin.MELLIN(npts=8,extended=True)
        conf['mellin-kaon']=mellin.MELLIN(npts=8,extended=True)
        conf['dmellin']= mellin.DMELLIN(nptsN=4,nptsM=4)
        conf['alphaS']= alphaS.ALPHAS()
        conf['eweak'] = eweak.EWEAK()
        conf['dsmf']  = deuterium.DEUTERON('%s/nuclib/grids/deuteron'%os.environ['FITPACK'])
        conf['hsmf']  = None#helium.HELIUM(root='%s/grids/grids-nuc-smearing/helium'%os.environ['FITPACK'])

        if 'pdf parametrization' in conf:

            if conf['pdf parametrization']==0: conf['pdf']= pdf0.PDF()
            if conf['pdf parametrization']==1: conf['pdf']= pdf1.PDF()
            if conf['pdf parametrization']==2: conf['pdf']= pdf2.PDF()
            if conf['pdf parametrization']==3: conf['pdf']= pdf3.PDF()
            if conf['pdf parametrization']==4: conf['pdf']= pdf4.PDF()
            if conf['pdf parametrization']==10:conf['pdf']= pdf10.PDF()

            conf['pdfA']  = conf['pdf']
            conf['pdfB']  = conf['pdf']

        if 'ht' in conf:
            if conf['ht']:
                conf['t3ppdf']= ht0.T3PPDF()
                #--0 for additive (proton/neutron), 1 for additive (proton/deuteron), 2 for multiplicative (proton/deuteron)
                if conf['ht parametrization']==0:
                  conf['t4F2']  = ht0.T4()
                  conf['t4FL']  = ht0.T4()
                  conf['t4F3']  = ht0.T4()
                  conf['t4g1']  = ht0.T4()
                if conf['ht parametrization']==1:
                  conf['t4F2']  = ht1.T4()
                  conf['t4FL']  = ht1.T4()
                  conf['t4F3']  = ht1.T4()
                  conf['t4g1']  = ht1.T4()
                if conf['ht parametrization']==2:
                  conf['t4F2']  = ht2.T4()
                  conf['t4FL']  = ht2.T4()
                  conf['t4F3']  = ht2.T4()
                  conf['t4g1']  = ht2.T4()

        if 'offshell' in conf:
            if conf['offshell']:
                #--0 for polynomial w/o sum rule, 1 for polynomial w/ sum rule, 2 for factorized polynomial w/ sum rule, 3 for factorized polynomial w/o sum rule
                if conf['offshell parametrization']==0:
                  conf['F2off']  = offshell0.OFFSHELL()
                  conf['g1off']  = offshell0.OFFSHELL()
                if conf['offshell parametrization']==1:
                  conf['F2off']  = offshell1.OFFSHELL()
                  conf['g1off']  = offshell1.OFFSHELL()
                if conf['offshell parametrization']==2:
                  conf['F2off']  = offshell2.OFFSHELL()
                  conf['g1off']  = offshell2.OFFSHELL()
                if conf['offshell parametrization']==3:
                  conf['F2off']  = offshell3.OFFSHELL() 
                  conf['g1off']  = offshell3.OFFSHELL() 
                if conf['offshell parametrization']==4:
                  conf['F2off']  = offshell4.OFFSHELL() 
                  conf['g1off']  = offshell4.OFFSHELL() 

        if 'ppdf parametrization' in conf:

            if conf['ppdf parametrization']==0: conf['ppdf']= ppdf0.PPDF()
            if conf['ppdf parametrization']==3: conf['ppdf']= ppdf3.PPDF()

            conf['ppdfA']  = conf['ppdf']
            conf['ppdfB']  = conf['ppdf']

        if 'ffpion parametrization' in conf:

            if conf['ffpion parametrization']==0: conf['ffpion']= ff0.FF()
            if conf['ffpion parametrization']==1: conf['ffpion']= ffpion1.FF()

        if 'ffkaon parametrization' in conf:

            if conf['ffkaon parametrization']==0: conf['ffkaon']= ff0.FF()
            if conf['ffkaon parametrization']==1: conf['ffkaon']= ffkaon1.FF()

        if 'ffhadron parametrization' in conf:
  
            if conf['ffhadron parametrization']=='0h': conf['ffhadron']= ff0h.FF()
            if conf['ffhadron parametrization']=='sum' and conf['ffpion parametrization']==0 and conf['ffkaon parametrization']==0: conf['ffhadron']= ffsum.FF() 
            if conf['ffhadron parametrization']=='sum2' and conf['ffpion parametrization']==0 and conf['ffkaon parametrization']==0: conf['ffhadron']= ffsum2.FF() 

        if 'pdf-pion parametrization' in conf:

            if conf['pdf-pion parametrization']==0: conf['pdf-pion']= pdfpion0.PDF()
            if conf['pdf-pion parametrization']==1: conf['pdf-pion']= pdfpion1.PDF()
            if conf['pdf-pion parametrization']=='sea': conf['pdf-pion']= pdfpionsea.PDF()

        if 'pdf-kaon parametrization' in conf:

            if conf['pdf-kaon parametrization']==0: conf['pdf-kaon']= pdfkaon.PDF()

    def setup_parallel(self,nworkers):
        self.parallel=PARALLEL()
        self.parallel.task=self.task
        self.parallel.set_state=self.set_state
        self.parallel.setup_master()
        self.parallel.setup_workers(nworkers)
        self.nworkers=nworkers
        self.requests=self.get_requests()

    def get_state(self):
        state={}
        if 'pdf'      in conf: state['pdf']      = conf['pdf'].get_state()
        if 'ppdf'     in conf: state['ppdf']     = conf['ppdf'].get_state()
        if 'ffpion'   in conf: state['ffpion']   = conf['ffpion'].get_state()
        if 'ffkaon'   in conf: state['ffkaon']   = conf['ffkaon'].get_state()
        if 'ffhadron' in conf: state['ffhadron'] = conf['ffhadron'].get_state()
        if 'pdf-pion' in conf: state['pdf-pion'] = conf['pdf-pion'].get_state()
        if 'ln'       in conf: state['ln']       = conf['ln'].get_state()
        if 'lh'       in conf: state['lh']       = conf['lh'].get_state()
        if 't4F2'     in conf: state['t4F2']     = conf['t4F2'].get_state()
        if 't4FL'     in conf: state['t4FL']     = conf['t4FL'].get_state()
        if 't4F3'     in conf: state['t4F3']     = conf['t4F3'].get_state()
        if 't4g1'     in conf: state['t4g1']     = conf['t4g1'].get_state()
        if 'F2off'    in conf: state['F2off']    = conf['F2off'].get_state()
        if 'g1off'    in conf: state['g1off']    = conf['g1off'].get_state()

        if 'zrap'     in conf: state['zrap']     = conf['zrap'].get_state()
        if 'wzrv'     in conf: state['wzrv']     = conf['wzrv'].get_state()
        return state

    def set_state(self,state):
        if 'pdf'       in conf: conf['pdf'].set_state(state['pdf'])
        if 'ppdf'      in conf: conf['ppdf'].set_state(state['ppdf'])
        if 'ffpion'    in conf: conf['ffpion'].set_state(state['ffpion'])
        if 'ffkaon'    in conf: conf['ffkaon'].set_state(state['ffkaon'])
        if 'ffhadron'  in conf: conf['ffhadron'].set_state(state['ffhadron'])
        if 'pdf-pion'  in conf: conf['pdf-pion'].set_state(state['pdf-pion'])
        if 'ln'        in conf: conf['ln'].set_state(state['ln'])
        if 'lh'        in conf: conf['lh'].set_state(state['lh'])
        if 't4F2'      in conf: conf['t4F2'].set_state(state['t4F2'])
        if 't4FL'      in conf: conf['t4FL'].set_state(state['t4FL'])
        if 't4F3'      in conf: conf['t4F3'].set_state(state['t4F3'])
        if 't4g1'      in conf: conf['t4g1'].set_state(state['t4g1'])
        if 'F2off'     in conf: conf['F2off'].set_state(state['F2off'])
        if 'g1off'     in conf: conf['g1off'].set_state(state['g1off'])

        if 'zrap'      in conf: conf['zrap'].set_state(state['zrap'])
        if 'wzrv'      in conf: conf['wzrv'].set_state(state['wzrv'])

    def distribute_requests(self,container,requests):
        cnt=0
        for request in requests:
            container[cnt].append(request)
            cnt+=1
            if cnt==self.nworkers: cnt=0

    def get_requests(self):
        container=[[] for _ in range(self.nworkers)]
        if 'idis'    in conf['datasets']:  self.distribute_requests(container,self.idisres.requests)
        if 'pidis'   in conf['datasets']:  self.distribute_requests(container,self.pidisres.requests)
        if 'sidis'   in conf['datasets']:  self.distribute_requests(container,self.sidisres.requests)
        if 'psidis'  in conf['datasets']:  self.distribute_requests(container,self.psidisres.requests)
        if 'dy'      in conf['datasets']:  self.distribute_requests(container,self.dyres.requests)
        if 'wzrv'    in conf['datasets']:  self.distribute_requests(container,self.wzrvres.requests)
        if 'wasym'   in conf['datasets']:  self.distribute_requests(container,self.wasymres.requests)
        if 'zrap'    in conf['datasets']:  self.distribute_requests(container,self.zrapres.requests)
        if 'sia'     in conf['datasets']:  self.distribute_requests(container,self.siares.requests)
        if 'jet'     in conf['datasets']:  self.distribute_requests(container,self.jetres.requests)
        if 'pjet'    in conf['datasets']:  self.distribute_requests(container,self.pjetres.requests)
        if 'qpdf'    in conf['datasets']:  self.distribute_requests(container,self.qpdfres.requests)
        if 'dy-pion' in conf['datasets']:  self.distribute_requests(container,self.dy_pion_res.requests)
        if 'ln'      in conf['datasets']:  self.distribute_requests(container,self.lnres.requests)
        if 'lh'      in conf['datasets']:  self.distribute_requests(container,self.lhres.requests)
        if 'pion_qT' in conf['datasets']:  self.distribute_requests(container,self.pion_qTres.requests)
        return container

    def task(self,request):
        for i in range(len(request)):
            if  request[i]['reaction']=='idis'    : self.idisres.process_request(request[i])
            if  request[i]['reaction']=='pidis'   : self.pidisres.process_request(request[i])
            if  request[i]['reaction']=='sidis'   : self.sidisres.process_request(request[i])
            if  request[i]['reaction']=='psidis'  : self.psidisres.process_request(request[i])
            if  request[i]['reaction']=='dy'      : self.dyres.process_request(request[i])
            if  request[i]['reaction']=='wzrv'    : self.wzrvres.process_request(request[i])
            if  request[i]['reaction']=='wasym'   : self.wasymres.process_request(request[i])
            if  request[i]['reaction']=='zrap'    : self.zrapres.process_request(request[i])
            if  request[i]['reaction']=='sia'     : self.siares.process_request(request[i])
            if  request[i]['reaction']=='jet'     : self.jetres.process_request(request[i])
            if  request[i]['reaction']=='pjet'    : self.pjetres.process_request(request[i])
            if  request[i]['reaction']=='qpdf'    : self.qpdfres.process_request(request[i])
            if  request[i]['reaction']=='dy-pion' : self.dy_pion_res.process_request(request[i])
            if  request[i]['reaction']=='ln'      : self.lnres.process_request(request[i])
            if  request[i]['reaction']=='lh'      : self.lhres.process_request(request[i])
            if  request[i]['reaction']=='pion_qT' : self.pion_qTres.process_request(request[i])
        return request

    def setup_idis(self):
        conf['idis tabs']=obslib.idis.reader.READER().load_data_sets('idis')
        conf['idis stfuncs']=obslib.idis.theory.STFUNCS()
        #if conf['hq']==True: conf['hqstfuncs']=obslib.dis.hqstfuncs.HQSTFUNCS()
        self.idisres=obslib.idis.residuals.RESIDUALS()

    def setup_sidis(self):
        conf['sidis tabs']=obslib.sidis.reader.READER().load_data_sets('sidis')
        conf['sidis stfuncs']=obslib.sidis.theory.STFUNCS()
        
        if 'idis stfuncs' not in conf:
            conf['idis stfuncs']=obslib.idis.theory.STFUNCS()
        
        self.sidisres=obslib.sidis.residuals.RESIDUALS()

    def setup_dy(self):
        conf['dy tabs']=obslib.dy.reader.READER().load_data_sets('dy')
        conf['dy']=obslib.dy.dymellin.DYMELLIN()
        conf['dy'].load_melltab()
        self.dyres=obslib.dy.residuals.RESIDUALS()

    def setup_wzrv(self):
        conf['wzrv tabs']=obslib.wzrv.reader.READER().load_data_sets('wzrv')
        conf['wzrv']=obslib.wzrv.wzmellin.WZMELLIN('dense',90)
        self.wzrvres=obslib.wzrv.residuals.RESIDUALS()

    def setup_wasym(self):
        conf['wasym tabs']=obslib.wasym.reader.READER().load_data_sets('wasym')
        conf['wasym']=obslib.wasym.wmellin.WMELLIN()
        conf['wasym'].load_melltab()
        self.wasymres=obslib.wasym.residuals.RESIDUALS()

    def setup_zrap(self):
        conf['zrap tabs']=obslib.zrap.reader.READER().load_data_sets('zrap')
        conf['zrap']=obslib.zrap.zmellin.ZMELLIN()
        self.zrapres=obslib.zrap.residuals.RESIDUALS()

    def setup_sia(self):
        conf['sia tabs']=obslib.sia.reader.READER().load_data_sets('sia')
        if 'ffpion parametrization'   in conf:  conf['sia pion']=obslib.sia.theory.SIA('pion')
        if 'ffkaon parametrization'   in conf:  conf['sia kaon']=obslib.sia.theory.SIA('kaon')
        if 'ffhadron parametrization' in conf:  conf['sia hadron']=obslib.sia.theory.SIA('hadron')
        self.siares=obslib.sia.residuals.RESIDUALS()

    def setup_pidis(self):
        conf['pidis tabs']=obslib.pidis.reader.READER().load_data_sets('pidis')
        conf['idis stfuncs']=obslib.idis.theory.STFUNCS()
        conf['idis pstfuncs']=obslib.pidis.theory.PSTFUNCS()
        self.pidisres=obslib.pidis.residuals.RESIDUALS()

    def setup_psidis(self):
        conf['psidis tabs']=obslib.psidis.reader.READER().load_data_sets('psidis')
        conf['sidis stfuncs']=obslib.sidis.theory.STFUNCS()
        conf['sidis pstfuncs']=obslib.psidis.theory.PSTFUNCS()
        self.psidisres=obslib.psidis.residuals.RESIDUALS()

    def setup_jet(self):
        if 'jet' not in conf:
            conf['jet'] = obslib.jets.jet_mellin.JETMELLIN()
        conf['jet tabs'] = obslib.jets.reader.READER().load_data_sets('jet')
        conf['jet'].load_factab(conf['jet_qr_fit']['method'], conf['jet_qr_fit']['f_scale'], conf['jet_qr_fit']['r_scale'])
        self.jetres = obslib.jets.residuals.RESIDUALS()

    def setup_pjet(self):
        conf['pjet'] = obslib.pjets.pjet_mellin.PJETMELLIN()
        conf['pjet tabs'] = obslib.pjets.reader.READER().load_data_sets('pjet')
        conf['pjet'].load_factab(conf['pjet_qr_fit']['method'], conf['pjet_qr_fit']['f_scale'], conf['pjet_qr_fit']['r_scale'])
        if 'jet' not in conf:
            conf['jet'] = obslib.jets.jet_mellin.JETMELLIN()
        conf['jet'].load_factab(conf['pjet_qr_fit']['method'], conf['pjet_qr_fit']['f_scale'], conf['pjet_qr_fit']['r_scale'], conf['pjet tabs'])
        self.pjetres = obslib.pjets.residuals.RESIDUALS()

    def setup_dy_pion(self):
        conf['path2dytab-hybrid']='%s/grids/grids-dypion/'%os.environ['FITPACK']
        conf['dy-pion tabs']=obslib.dy.reader.READER_PIONS().load_data_sets('dy-pion')
        conf['dy-pion']=obslib.dy.piontheory.DY_PION()
        conf['dy-pion'].load_melltab_hybrid()
        self.dy_pion_res=obslib.dy.residuals.RESIDUALS_PIONS()

    def setup_ln(self):
        conf['ln tabs']=obslib.ln.reader.READER().load_data_sets('ln')
        conf['pion-stfuncs']=obslib.ln.theory.STFUNCS()
        conf['ln']=obslib.ln.theory.LN()
        if 'idis stfuncs' not in conf:
            conf['path2idistab']='%s/grids/grids-idis/distab'%(os.environ['FITPACK'])
            conf['idis stfuncs']=obslib.idis.theory.STFUNCS()
        self.lnres=obslib.ln.residuals.RESIDUALS()

    def setup_lh(self):
        conf['lh tabs']=obslib.lh.reader.READER().load_data_sets('lh')
        conf['kaon-stfuncs']=obslib.lh.theory.STFUNCS()
        conf['lh']=obslib.lh.theory.LH()
        if 'idis stfuncs' not in conf:
            conf['path2idistab']='%s/grids/grids-idis/distab'%(os.environ['FITPACK'])
            conf['idis stfuncs']=obslib.idis.theory.STFUNCS()
        self.lhres=obslib.lh.residuals.RESIDUALS()

    def setup_pion_qT(self):
        conf['pion_qT tabs']=obslib.pion_qT.reader.READER().load_data_sets('pion_qT')
        #conf['pion_qT']=obslib.pion_qT.new_y_theory.PION_QT()
        conf['pion_qT']=obslib.pion_qT.theory.PION_QT()
        conf['pion_qT-asy']=obslib.pion_qT.theoryasy.PION_QT()
        self.pion_qTres=obslib.pion_qT.residuals.RESIDUALS()

    def setup_qpdf(self):
        conf['qpdf tabs']=obslib.qpdf.reader.READER().load_data_sets('qpdf')
        conf['qpdf']=obslib.qpdf.theory.QPDF()
        conf['qpdf'].load_mdata()
        self.qpdfres=obslib.qpdf.residuals.RESIDUALS()

    def get_residuals(self,par):
        self.parman.set_new_params(par)
        state=self.get_state()
        self.parallel.update_workers(state)
        results=self.parallel.send_tasks(self.requests)

        #--update tables with the new theory values
        for chunk in results:
            for request in chunk:
                if request['reaction']=='idis'   : self.idisres.update_tabs_external(request)
                if request['reaction']=='pidis'  : self.pidisres.update_tabs_external(request)
                if request['reaction']=='sidis'  : self.sidisres.update_tabs_external(request)
                if request['reaction']=='psidis' : self.psidisres.update_tabs_external(request)
                if request['reaction']=='dy'     : self.dyres.update_tabs_external(request)
                if request['reaction']=='wzrv'   : self.wzrvres.update_tabs_external(request)
                if request['reaction']=='wasym'  : self.wasymres.update_tabs_external(request)
                if request['reaction']=='zrap'   : self.zrapres.update_tabs_external(request)
                if request['reaction']=='sia'    : self.siares.update_tabs_external(request)
                if request['reaction']=='jet'    : self.jetres.update_tabs_external(request)
                if request['reaction']=='pjet'   : self.pjetres.update_tabs_external(request)
                if request['reaction']=='qpdf'   : self.qpdfres.update_tabs_external(request)
                if request['reaction']=='dy-pion': self.dy_pion_res.update_tabs_external(request)
                if request['reaction']=='ln'     : self.lnres.update_tabs_external(request)
                if request['reaction']=='lh'     : self.lhres.update_tabs_external(request)
                if request['reaction']=='pion_qT': self.pion_qTres.update_tabs_external(request)

        if self.parallel_mode: calc=False
        else: calc=True

        #--compute residuals
        res,rres,nres=[],[],[]
        if 'idis' in conf['datasets']:
            out=self.idisres.get_residuals(calc=calc)
            res=np.append(res,out[0])
            rres=np.append(rres,out[1])
            nres=np.append(nres,out[2])
        if 'pidis' in conf['datasets']:
            out=self.pidisres.get_residuals(calc=calc)
            res=np.append(res,out[0])
            rres=np.append(rres,out[1])
            nres=np.append(nres,out[2])
        if 'sidis' in conf['datasets']:
            out=self.sidisres.get_residuals(calc=calc)
            res=np.append(res,out[0])
            rres=np.append(rres,out[1])
            nres=np.append(nres,out[2])
        if 'psidis' in conf['datasets']:
            out=self.psidisres.get_residuals(calc=calc)
            res=np.append(res,out[0])
            rres=np.append(rres,out[1])
            nres=np.append(nres,out[2])
        if 'dy'   in conf['datasets']:
            out=self.dyres.get_residuals(calc=calc)
            res=np.append(res,out[0])
            rres=np.append(rres,out[1])
            nres=np.append(nres,out[2])
        if 'wzrv'   in conf['datasets']:
            out=self.wzrvres.get_residuals(calc=calc)
            res=np.append(res,out[0])
            rres=np.append(rres,out[1])
            nres=np.append(nres,out[2])
        if 'wasym'   in conf['datasets']:
            out=self.wasymres.get_residuals(calc=calc)
            res=np.append(res,out[0])
            rres=np.append(rres,out[1])
            nres=np.append(nres,out[2])
        if 'zrap'   in conf['datasets']:
            out=self.zrapres.get_residuals(calc=calc)
            res=np.append(res,out[0])
            rres=np.append(rres,out[1])
            nres=np.append(nres,out[2])
        if 'sia' in conf['datasets']:
            out=self.siares.get_residuals(calc=calc)
            res=np.append(res,out[0])
            rres=np.append(rres,out[1])
            nres=np.append(nres,out[2])
        if 'jet' in conf['datasets']:
            out=self.jetres.get_residuals(calc=calc)
            res=np.append(res,out[0])
            rres=np.append(rres,out[1])
            nres=np.append(nres,out[2])
        if 'pjet' in conf['datasets']:
            out=self.pjetres.get_residuals(calc=calc)
            res=np.append(res,out[0])
            rres=np.append(rres,out[1])
            nres=np.append(nres,out[2])
        if 'qpdf' in conf['datasets']:
            out=self.qpdfres.get_residuals(calc=calc)
            res=np.append(res,out[0])
            rres=np.append(rres,out[1])
            nres=np.append(nres,out[2])
        if 'dy-pion' in conf['datasets']:
            out=self.dy_pion_res.get_residuals(calc=calc)
            res=np.append(res,out[0])
            rres=np.append(rres,out[1])
            nres=np.append(nres,out[2])
            if 'penalty' in conf and conf['penalty']==True: res=np.append(res,parman.check_residual_penalty())
        if 'ln' in conf['datasets']:
            out=self.lnres.get_residuals(calc=calc)
            res=np.append(res,out[0])
            rres=np.append(rres,out[1])
            nres=np.append(nres,out[2])
        if 'lh' in conf['datasets']:
            out=self.lhres.get_residuals(calc=calc)
            res=np.append(res,out[0])
            rres=np.append(rres,out[1])
            nres=np.append(nres,out[2])
        if 'pion_qT' in conf['datasets']:
            out=self.pion_qTres.get_residuals(calc=calc)
            res=np.append(res,out[0])
            rres=np.append(rres,out[1])
            nres=np.append(nres,out[2])
        return res,rres,nres

    def get_data_info(self):

        #--compute residuals
        reaction=[]
        if 'idis' in conf['datasets']:
            out=self.idisres.get_residuals(calc=False)
            reaction.extend(['idis' for _ in out[0]])
        if 'pidis' in conf['datasets']:
            out=self.pidisres.get_residuals(calc=False)
            reaction.extend(['pidis' for _ in out[0]])
        if 'sidis' in conf['datasets']:
            out=self.sidisres.get_residuals(calc=False)
            reaction.extend(['sidis' for _ in out[0]])
        if 'psidis' in conf['datasets']:
            out=self.psidisres.get_residuals(calc=False)
            reaction.extend(['psidis' for _ in out[0]])
        if 'dy'   in conf['datasets']:
            out=self.dyres.get_residuals(calc=False)
            reaction.extend(['dy' for _ in out[0]])
        if 'wzrv'   in conf['datasets']:
            out=self.wzrvres.get_residuals(calc=False)
            reaction.extend(['wzrv' for _ in out[0]])
        if 'wasym'   in conf['datasets']:
            out=self.wasymres.get_residuals(calc=False)
            reaction.extend(['wasym' for _ in out[0]])
        if 'zrap'   in conf['datasets']:
            out=self.zrapres.get_residuals(calc=False)
            reaction.extend(['zrap' for _ in out[0]])
        if 'sia' in conf['datasets']:
            out=self.siares.get_residuals(calc=False)
            reaction.extend(['sia' for _ in out[0]])
        if 'jet' in conf['datasets']:
            out=self.jetres.get_residuals(calc=False)
            reaction.extend(['jet' for _ in out[0]])
        if 'pjet' in conf['datasets']:
            out=self.pjetres.get_residuals(calc=False)
            reaction.extend(['pjet' for _ in out[0]])
        if 'qpdf' in conf['datasets']:
            out=self.qpdfres.get_residuals(calc=False)
            reaction.extend(['qpdf' for _ in out[0]])
        if 'dy-pion' in conf['datasets']:
            out=self.dy_pion_res.get_residuals(calc=False)
            reaction.extend(['dy-pion' for _ in out[0]])
        if 'ln'  in conf['datasets']:
            out=self.lnres.get_residuals(calc=False)
            reaction.extend(['ln' for _ in out[0]])
        if 'lh'  in conf['datasets']:
            out=self.lhres.get_residuals(calc=False)
            reaction.extend(['lh' for _ in out[0]])
        if 'pion_qT'  in conf['datasets']:
            out=self.pion_qTres.get_residuals(calc=False)
            reaction.extend(['pion_qT' for _ in out[0]])
        return reaction

    def gen_report(self,verb=0,level=0):
        L=[]
        if 'idis'    in conf['datasets']: L.extend(self.idisres.gen_report(verb,level))
        if 'pidis'   in conf['datasets']: L.extend(self.pidisres.gen_report(verb,level))
        if 'sidis'   in conf['datasets']: L.extend(self.sidisres.gen_report(verb,level))
        if 'psidis'  in conf['datasets']: L.extend(self.psidisres.gen_report(verb,level))
        if 'dy'      in conf['datasets']: L.extend(self.dyres.gen_report(verb,level))
        if 'wzrv'    in conf['datasets']: L.extend(self.wzrvres.gen_report(verb,level))
        if 'wasym'   in conf['datasets']: L.extend(self.wasymres.gen_report(verb,level))
        if 'zrap'    in conf['datasets']: L.extend(self.zrapres.gen_report(verb,level))
        if 'sia'     in conf['datasets']: L.extend(self.siares.gen_report(verb,level))
        if 'jet'     in conf['datasets']: L.extend(self.jetres.gen_report(verb,level))
        if 'pjet'    in conf['datasets']: L.extend(self.pjetres.gen_report(verb,level))
        if 'qpdf'    in conf['datasets']: L.extend(self.qpdfres.gen_report(verb,level))
        if 'dy-pion' in conf['datasets']: L.extend(self.dy_pion_res.gen_report(verb,level))
        if 'ln'      in conf['datasets']: L.extend(self.lnres.gen_report(verb,level))
        if 'lh'      in conf['datasets']: L.extend(self.lhres.gen_report(verb,level))
        if 'pion_qT' in conf['datasets']: L.extend(self.pion_qTres.gen_report(verb,level))
        return L

    def get_chi2(self):

        data={}
        if 'idis'    in conf['datasets']: data.update(self.idisres.get_chi2())
        if 'pidis'   in conf['datasets']: data.update(self.pidisres.get_chi2())
        if 'sidis'   in conf['datasets']: data.update(self.sidisres.get_chi2())
        if 'psidis'  in conf['datasets']: data.update(self.psidisres.get_chi2())
        if 'dy'      in conf['datasets']: data.update(self.dyres.get_chi2())
        if 'wzrv'    in conf['datasets']: data.update(self.wzrvres.get_chi2())
        if 'wasym'   in conf['datasets']: data.update(self.wasymres.get_chi2())
        if 'zrap'    in conf['datasets']: data.update(self.zrapres.get_chi2())
        if 'sia'     in conf['datasets']: data.update(self.siares.get_chi2())
        if 'jet'     in conf['datasets']: data.update(self.jetres.get_chi2())
        if 'pjet'    in conf['datasets']: data.update(self.pjetres.get_chi2())
        if 'qpdf'    in conf['datasets']: data.update(self.qpdfres.get_chi2())
        if 'dy-pion' in conf['datasets']: data.update(self.dy_pion_res.get_chi2())
        if 'ln'      in conf['datasets']: data.update(self.lnres.get_chi2())
        if 'lh'      in conf['datasets']: data.update(self.lhres.get_chi2())
        if 'pion_qT' in conf['datasets']: data.update(self.pion_qTres.get_chi2())
        return data

    def test(self,ntasks=1):
        #--loop over states
        print('='*20)
        t=time.time()
        for _ in range(ntasks):
            par=self.parman.gen_flat(setup=True)
            res,rres,nres=self.get_residuals(par)
            chi2=np.sum(res**2)
            print '(%d/%d) chi2=%f'%(_ + 1, ntasks, chi2)
        print('='*20)
        elapsed_time=time.time()-t
        print('with %d workers, elapsed time: %f' % (self.nworkers, elapsed_time))
        return elapsed_time

    def shutdown(self):
        if self.parallel_mode:
            self.parallel.stop_workers()


