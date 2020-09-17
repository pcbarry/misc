import os
conf={}

#--setup posterior sampling

conf['bootstrap']=False
conf['flat par']=False
conf['ftol']=1e-8

#--setup qcd evolution

conf['dglap mode']='truncated'
conf['alphaSmode']='backward'
conf['order'] = 'NLO'
conf['Q20']   = 1.27**2

#--setup for idis

conf['tmc']   = False
conf['ht']    = False
conf['nuc smearing']=True
conf['sidis nuc smearing']=False
conf['hq']=False

#--grids

conf['path2idistab']   = '%s/grids/grids-idis/distab/'%os.environ['FITPACK']
conf['path2pidistab']   = '%s/grids/grids-pidis/distab/'%os.environ['FITPACK']

#--datasets

conf['datasets']={}

#--lepton-hadron reactions

Q2cut=1.3**2
W2cut=10.0

##--IDIS
conf['datasets']['idis']={}
conf['datasets']['idis']['filters']=[]
conf['datasets']['idis']['filters'].append("Q2>%f"%Q2cut)
conf['datasets']['idis']['filters'].append("W2>%f"%W2cut)
conf['datasets']['idis']['xlsx']={}
conf['datasets']['idis']['xlsx'][10010]='idis/expdata/10010.xlsx' # proton   | F2            | SLAC
conf['datasets']['idis']['xlsx'][10011]='idis/expdata/10011.xlsx' # deuteron | F2            | SLAC
conf['datasets']['idis']['xlsx'][10016]='idis/expdata/10016.xlsx' # proton   | F2            | BCDMS
conf['datasets']['idis']['xlsx'][10017]='idis/expdata/10017.xlsx' # deuteron | F2            | BCDMS
conf['datasets']['idis']['xlsx'][10020]='idis/expdata/10020.xlsx' # proton   | F2            | NMC
conf['datasets']['idis']['xlsx'][10021]='idis/expdata/10021.xlsx' # d/p      | F2d/F2p       | NMC
conf['datasets']['idis']['norm']={}
conf['datasets']['idis']['norm'][10010]={'value':   1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['idis']['norm'][10011]={'value':   1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['idis']['norm'][10016]={'value':   1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['idis']['norm'][10017]={'value':   1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['idis']['norm'][10020]={'value':   1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}


#--parameters
conf['params'] = {}

#--pdf parameters

conf['pdf parametrization'] = 2
conf['params']['pdf'] = {}

conf['params']['pdf']['g1 N']    ={'value':    1, 'min':  None, 'max':  None, 'fixed': True }
conf['params']['pdf']['g1 a']    ={'value': -0.5, 'min':  -1.9, 'max':     1, 'fixed': False}
conf['params']['pdf']['g1 b']    ={'value':    6, 'min':     0, 'max':    10, 'fixed': False}

conf['params']['pdf']['uv1 N']   ={'value':    1, 'min':  None, 'max':  None, 'fixed': True }
conf['params']['pdf']['uv1 a']   ={'value': -0.5, 'min':  -0.5, 'max':     1, 'fixed': False}
conf['params']['pdf']['uv1 b']   ={'value':    6, 'min':     0, 'max':    10, 'fixed': False}

conf['params']['pdf']['dv1 N']   ={'value':    1, 'min':  None, 'max':  None, 'fixed': True }
conf['params']['pdf']['dv1 a']   ={'value': -0.5, 'min':  -0.5, 'max':     1, 'fixed': False}
conf['params']['pdf']['dv1 b']   ={'value':    6, 'min':     0, 'max':    10, 'fixed': False}

conf['params']['pdf']['db1 N']   ={'value':    1, 'min':     0, 'max':     1, 'fixed': False}
conf['params']['pdf']['db1 a']   ={'value': -0.5, 'min':    -1, 'max':     1, 'fixed': False}
conf['params']['pdf']['db1 b']   ={'value':    6, 'min':     0, 'max':    10, 'fixed': False}

conf['params']['pdf']['ub1 N']   ={'value':    1, 'min':     0, 'max':     1, 'fixed': False}
conf['params']['pdf']['ub1 a']   ={'value': -0.5, 'min':    -1, 'max':     1, 'fixed': False}
conf['params']['pdf']['ub1 b']   ={'value':    6, 'min':     0, 'max':    10, 'fixed': False}

conf['params']['pdf']['s1 N']    ={'value':    1, 'min':     0, 'max':     1, 'fixed': False}
conf['params']['pdf']['s1 a']    ={'value':    0, 'min':    -1, 'max':     1, 'fixed': False}
conf['params']['pdf']['s1 b']    ={'value':    6, 'min':     0, 'max':    10, 'fixed': False}

conf['params']['pdf']['sb1 N']   ={'value':    1, 'min':     0, 'max':     1, 'fixed': 's1 N'}
conf['params']['pdf']['sb1 a']   ={'value':    0, 'min':    -1, 'max':     1, 'fixed': 's1 a'}
conf['params']['pdf']['sb1 b']   ={'value':    6, 'min':     0, 'max':    10, 'fixed': 's1 b'}

conf['params']['pdf']['sea1 N']  ={'value':  0.5, 'min':     0, 'max':     1, 'fixed': False}
conf['params']['pdf']['sea1 a']  ={'value': -1.5, 'min':  -1.9, 'max':    -1, 'fixed': False}
conf['params']['pdf']['sea1 b']  ={'value':    6, 'min':     0, 'max':    10, 'fixed': False}

conf['params']['pdf']['sea2 N']  ={'value':    1, 'min':     0, 'max':     1, 'fixed': 'sea1 N'}
conf['params']['pdf']['sea2 a']  ={'value': -1.5, 'min':  -1.9, 'max':    -1, 'fixed': 'sea1 a'}
conf['params']['pdf']['sea2 b']  ={'value':    6, 'min':     0, 'max':    10, 'fixed': 'sea1 b'}

#--steps
conf['steps']={}

#--idis and no hera
conf['steps'][1]={}
conf['steps'][1]['dep']=[]
conf['steps'][1]['active distributions']=['pdf']
conf['steps'][1]['datasets']={}
conf['steps'][1]['datasets']['idis']=[]
conf['steps'][1]['datasets']['idis'].append(10010) # proton   | F2            | SLAC
conf['steps'][1]['datasets']['idis'].append(10011) # deuteron | F2            | SLAC
conf['steps'][1]['datasets']['idis'].append(10016) # proton   | F2            | BCDMS
conf['steps'][1]['datasets']['idis'].append(10017) # deuteron | F2            | BCDMS
conf['steps'][1]['datasets']['idis'].append(10020) # proton   | F2            | NMC
conf['steps'][1]['datasets']['idis'].append(10021) # d/p      | F2d/F2p       | NMC




