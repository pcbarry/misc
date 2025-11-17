import os

conf={}
conf['ftol']=1e-6
conf['bootstrap']=False
conf['flat par']=True

# Data
conf['datasets']={}
conf['datasets']['dvcs']={}
conf['datasets']['dvcs']['filters']=[]
conf['datasets']['dvcs']['filters'].append('Q2>1.7')
conf['datasets']['dvcs']['filters'].append('-t<0.2*Q2')
conf['datasets']['dvcs']['xlsx']={}
# conf['datasets']['dvcs']['xlsx'][10001]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_1_Table_1.xlsx' # p ALU HERMES 2001
# conf['datasets']['dvcs']['xlsx'][10002]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_2_Table_1.xlsx' # p A_LU_sinphi CLAS 2001
# conf['datasets']['dvcs']['xlsx'][10003]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_2_Table_2.xlsx' # p A_LU_sin2phi CLAS 2001
# conf['datasets']['dvcs']['xlsx'][10004]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_3_Table_1a.xlsx' # p dSig/dt H1 1996-1997
# conf['datasets']['dvcs']['xlsx'][10005]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_3_Table_1b.xlsx' # p dSig/dt H1 1999-2000
# conf['datasets']['dvcs']['xlsx'][10006]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_4_Table_1.xlsx' # p A_UL_sinphi CLAS 2006
# conf['datasets']['dvcs']['xlsx'][10007]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_4_Table_2.xlsx' # p A_UL_sin2phi CLAS 2006
# conf['datasets']['dvcs']['xlsx'][10008]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_6_Table_1a.xlsx' # p A_C_cosphi HERMES 2007
# conf['datasets']['dvcs']['xlsx'][10009]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_7_Table_1.xlsx' # p ALU CLAS 2008
# conf['datasets']['dvcs']['xlsx'][10010]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_8_Table_1axB.xlsx' # p A_C_cos0phi HERMES 2008
# conf['datasets']['dvcs']['xlsx'][10011]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_8_Table_1bxB.xlsx' # p A_C_cosphi HERMES 2008
# conf['datasets']['dvcs']['xlsx'][10012]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_8_Table_1cxB.xlsx' # p A_UT_DVCS_cos0phi HERMES 2008
# conf['datasets']['dvcs']['xlsx'][10013]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_8_Table_1dxB.xlsx' # p A_UT_I_cos0phi HERMES 2008
# conf['datasets']['dvcs']['xlsx'][10014]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_8_Table_1exB.xlsx' # p A_UT_I_cosphi HERMES 2008
# conf['datasets']['dvcs']['xlsx'][10015]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_8_Table_1fxB.xlsx' # p A_UT_I_sinphi HERMES 2008
# conf['datasets']['dvcs']['xlsx'][10028]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_9_Table_4.xlsx' # p dSig/dt ZEUS 2009
# conf['datasets']['dvcs']['xlsx'][10029]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_10_Table_1.xlsx' # p ALU CLAS 2009
# conf['datasets']['dvcs']['xlsx'][10030]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_10_Table_2.xlsx' # p ALU CLAS 2009
# conf['datasets']['dvcs']['xlsx'][10031]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_11_Table_3a.xlsx' # p dSig/dt H1 2009
# conf['datasets']['dvcs']['xlsx'][10032]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_11_Table_3b.xlsx' # p dSig/dt H1 2009
# conf['datasets']['dvcs']['xlsx'][10034]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_12_Table_2_2.xlsx' # p A_LU_I_sinphi HERMES 2009
# conf['datasets']['dvcs']['xlsx'][10035]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_12_Table_2_3.xlsx' # p A_LU_I_sin2phi HERMES 2009
# conf['datasets']['dvcs']['xlsx'][10036]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_12_Table_4_0.xlsx' # p A_C_cos0phi HERMES 2009
# conf['datasets']['dvcs']['xlsx'][10037]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_12_Table_4_1.xlsx' # p A_C_cosphi HERMES 2009
# conf['datasets']['dvcs']['xlsx'][10038]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_12_Table_4_2.xlsx' # p A_C_cos2phi HERMES 2009
# conf['datasets']['dvcs']['xlsx'][10039]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_12_Table_4_3.xlsx' # p A_C_cos3phi HERMES 2009
# conf['datasets']['dvcs']['xlsx'][10052]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_13_Table_4axB.xlsx' # p A_UL_sinphi HERMES 2010
# conf['datasets']['dvcs']['xlsx'][10053]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_13_Table_4bxB.xlsx' # p A_UL_sin2phi HERMES 2010
# conf['datasets']['dvcs']['xlsx'][10054]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_13_Table_4cxB.xlsx' # p A_UL_sin3phi HERMES 2010
# conf['datasets']['dvcs']['xlsx'][10055]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_13_Table_4dxB.xlsx' # p A_LL_cos0phi HERMES 2010
# conf['datasets']['dvcs']['xlsx'][10056]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_13_Table_4exB.xlsx' # p A_LL_cosphi HERMES 2010
# conf['datasets']['dvcs']['xlsx'][10057]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_13_Table_4fxB.xlsx' # p A_LL_cos2phi HERMES 2010
# conf['datasets']['dvcs']['xlsx'][10058]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_14_Table_2a.xlsx' # p A_LT_I_cos0phi HERMES 2011
# conf['datasets']['dvcs']['xlsx'][10059]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_14_Table_2b.xlsx' # p A_LT_I_cosphi HERMES 2011
# conf['datasets']['dvcs']['xlsx'][10060]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_14_Table_2c.xlsx' # p A_LT_I_cos2phi HERMES 2011
# conf['datasets']['dvcs']['xlsx'][10061]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_14_Table_2d.xlsx' # p A_LT_I_sinphi HERMES 2011
# conf['datasets']['dvcs']['xlsx'][10062]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_14_Table_2e.xlsx' # p A_LT_I_sin2phi HERMES 2011
# conf['datasets']['dvcs']['xlsx'][10063]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_14_Table_3a.xlsx' # p A_LT_DVCS_cos0phi HERMES 2011
# conf['datasets']['dvcs']['xlsx'][10064]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_14_Table_3b.xlsx' # p A_LT_DVCS_cosphi HERMES 2011
# conf['datasets']['dvcs']['xlsx'][10065]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_14_Table_3c.xlsx' # p A_LT_DVCS_sinphi HERMES 2011
# conf['datasets']['dvcs']['xlsx'][10083]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_15_Table_1b.xlsx' # p A_LU_I_sinphi HERMES 2012
# conf['datasets']['dvcs']['xlsx'][10084]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_15_Table_1c.xlsx' # p A_LU_I_sin2phi HERMES 2012
# conf['datasets']['dvcs']['xlsx'][10085]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_15_Table_2a.xlsx' # p A_C_cos0phi HERMES 2012
# conf['datasets']['dvcs']['xlsx'][10086]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_15_Table_2b.xlsx' # p A_C_cosphi HERMES 2012
# conf['datasets']['dvcs']['xlsx'][10087]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_15_Table_2c.xlsx' # p A_C_cos2phi HERMES 2012
# conf['datasets']['dvcs']['xlsx'][10088]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_15_Table_2d.xlsx' # p A_C_cos3phi HERMES 2012
# conf['datasets']['dvcs']['xlsx'][10103]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_16_Table_1.xlsx' # p ALU CLAS 2015
# conf['datasets']['dvcs']['xlsx'][10104]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_16_Table_2.xlsx' # p AUL CLAS 2015
# conf['datasets']['dvcs']['xlsx'][10105]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_16_Table_3.xlsx' # p ALL CLAS 2015
conf['datasets']['dvcs']['xlsx'][10106]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_17_UU.xlsx' # p d4SigUU Hall A 2015
conf['datasets']['dvcs']['xlsx'][10107]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_17_LU.xlsx' # p deltad4SigLU Hall A 2015
conf['datasets']['dvcs']['xlsx'][10108]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_18_Table_1.xlsx' # p d4SigUU CLAS 2015
conf['datasets']['dvcs']['xlsx'][10109]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_18_Table_2.xlsx' # p deltad4SigLU CLAS 2015
conf['datasets']['dvcs']['xlsx'][10110]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_19_Table_1.xlsx' # p deltad4SigLU Hall A 2017
conf['datasets']['dvcs']['xlsx'][10111]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_19_Table_2.xlsx' # p deltad4SigLU Hall A 2017
# conf['datasets']['dvcs']['xlsx'][10112]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_19_Table_3.xlsx' # p deltad4SigLU Hall A 2017
conf['datasets']['dvcs']['xlsx'][10113]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_19_Table_4.xlsx' # p d4SigUU Hall A 2017
# conf['datasets']['dvcs']['xlsx'][10114]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_19_Table_5.xlsx' # p d4SigUU Hall A 2017
conf['datasets']['dvcs']['xlsx'][10115]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_19_Table_6.xlsx' # p d4SigUU Hall A 2017
# conf['datasets']['dvcs']['xlsx'][10116]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_21_Table_1.xlsx' # p dSig/dt COMPASS 2019
# conf['datasets']['dvcs']['xlsx'][10117]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_24_UU.xlsx' # p d4SigUU Hall A 2022
# conf['datasets']['dvcs']['xlsx'][10118]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_24_LU.xlsx' # p deltad4SigLU Hall A 2022
# conf['datasets']['dvcs']['xlsx'][10119]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_25_01.xlsx' # p ALU CLAS 2023
# conf['datasets']['dvcs']['xlsx'][10120]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_25_02.xlsx' # p ALU CLAS 2023
conf['datasets']['dvcs']['xlsx'][10121]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_28_Table_1.xlsx' # p d4SigUU CLAS 2018
conf['datasets']['dvcs']['xlsx'][10122]='/work/JAM/moffatea/QGTDatabase/Tables/Ref_28_Table_2.xlsx' # p deltad4SigLU CLAS 2018
conf['datasets']['dvcs']['norm']={}
# conf['datasets']['dvcs']['norm'][10009]={'value': 1.0,'min': 0.8,'max': 1.2,'fixed': False,'forcemin':True,'forcemax':True}
# conf['datasets']['dvcs']['norm'][10012]={'value': 1.0,'min': 0.8,'max': 1.2,'fixed': False,'forcemin':True,'forcemax':True}
# conf['datasets']['dvcs']['norm'][10013]={'value': 1.0,'min': 0.8,'max': 1.2,'fixed': False,'forcemin':True,'forcemax':True}
# conf['datasets']['dvcs']['norm'][10014]={'value': 1.0,'min': 0.8,'max': 1.2,'fixed': False,'forcemin':True,'forcemax':True}
# conf['datasets']['dvcs']['norm'][10015]={'value': 1.0,'min': 0.8,'max': 1.2,'fixed': False,'forcemin':True,'forcemax':True}
# conf['datasets']['dvcs']['norm'][10028]={'value': 1.0,'min': 0.8,'max': 1.2,'fixed': False,'forcemin':True,'forcemax':True}
# conf['datasets']['dvcs']['norm'][10029]={'value': 1.0,'min': 0.8,'max': 1.2,'fixed': False,'forcemin':True,'forcemax':True}
# conf['datasets']['dvcs']['norm'][10030]={'value': 1.0,'min': 0.8,'max': 1.2,'fixed': False,'forcemin':True,'forcemax':True}
# conf['datasets']['dvcs']['norm'][10052]={'value': 1.0,'min': 0.8,'max': 1.2,'fixed': False,'forcemin':True,'forcemax':True}
# conf['datasets']['dvcs']['norm'][10053]={'value': 1.0,'min': 0.8,'max': 1.2,'fixed': False,'forcemin':True,'forcemax':True}
# conf['datasets']['dvcs']['norm'][10054]={'value': 1.0,'min': 0.8,'max': 1.2,'fixed': False,'forcemin':True,'forcemax':True}
# conf['datasets']['dvcs']['norm'][10055]={'value': 1.0,'min': 0.8,'max': 1.2,'fixed': False,'forcemin':True,'forcemax':True}
# conf['datasets']['dvcs']['norm'][10056]={'value': 1.0,'min': 0.8,'max': 1.2,'fixed': False,'forcemin':True,'forcemax':True}
# conf['datasets']['dvcs']['norm'][10057]={'value': 1.0,'min': 0.8,'max': 1.2,'fixed': False,'forcemin':True,'forcemax':True}
# conf['datasets']['dvcs']['norm'][10083]={'value': 1.0,'min': 0.8,'max': 1.2,'fixed': False,'forcemin':True,'forcemax':True}
# conf['datasets']['dvcs']['norm'][10084]={'value': 1.0,'min': 0.8,'max': 1.2,'fixed': False,'forcemin':True,'forcemax':True}
conf['datasets']['dvcs']['norm'][10106]={'value': 1.0,'min': 0.8,'max': 1.2,'fixed': False,'forcemin':True,'forcemax':True}
conf['datasets']['dvcs']['norm'][10107]={'value': 1.0,'min': 0.8,'max': 1.2,'fixed': False,'forcemin':True,'forcemax':True}
conf['datasets']['dvcs']['norm'][10110]={'value': 1.0,'min': 0.8,'max': 1.2,'fixed': False,'forcemin':True,'forcemax':True}
conf['datasets']['dvcs']['norm'][10111]={'value': 1.0,'min': 0.8,'max': 1.2,'fixed': False,'forcemin':True,'forcemax':True}
conf['datasets']['dvcs']['norm'][10113]={'value': 1.0,'min': 0.8,'max': 1.2,'fixed': False,'forcemin':True,'forcemax':True}
conf['datasets']['dvcs']['norm'][10115]={'value': 1.0,'min': 0.8,'max': 1.2,'fixed': False,'forcemin':True,'forcemax':True}

# Parameters
conf['params']={}

# Unpolarized PDFs
conf['params']['pdf']={}
conf['params']['pdf']['g1 N']={'value': 0.388,'fixed':True}
conf['params']['pdf']['g1 a']={'value': -0.623,'fixed':True}
conf['params']['pdf']['g1 b']={'value': 9.26,'fixed':True}
conf['params']['pdf']['g1 c']={'value': 0.0,'fixed':True}
conf['params']['pdf']['g1 d']={'value': 0.0,'fixed':True}
    
conf['params']['pdf']['uv1 N']={'value': 0.346,'fixed':True}
conf['params']['pdf']['uv1 a']={'value': -0.122,'fixed':True}
conf['params']['pdf']['uv1 b']={'value': 3.21,'fixed':True}
conf['params']['pdf']['uv1 c']={'value': 0.0,'fixed':True}
conf['params']['pdf']['uv1 d']={'value': 0.0,'fixed':True}
        
conf['params']['pdf']['dv1 N']={'value': 0.152,'fixed':True}
conf['params']['pdf']['dv1 a']={'value': -0.24,'fixed':True}
conf['params']['pdf']['dv1 b']={'value': 3.84,'fixed':True}
conf['params']['pdf']['dv1 c']={'value': 0.0,'fixed':True}
conf['params']['pdf']['dv1 d']={'value': 0.0,'fixed':True}
        
conf['params']['pdf']['mix N']={'value': 0.0,'fixed':True}
conf['params']['pdf']['mix a']={'value': 2.0,'fixed':True}
        
conf['params']['pdf']['db1 N']={'value': 3.68,'fixed':True}
conf['params']['pdf']['db1 a']={'value': -8.41,'fixed':True}
conf['params']['pdf']['db1 b']={'value': 5.31,'fixed':True}
conf['params']['pdf']['db1 c']={'value': 0.0,'fixed':True}
conf['params']['pdf']['db1 d']={'value': 0.0,'fixed':True}
    
conf['params']['pdf']['ub1 N']={'value': 1.95,'fixed':True}
conf['params']['pdf']['ub1 a']={'value': -9.94,'fixed':True}
conf['params']['pdf']['ub1 b']={'value': 8.39,'fixed':True}
conf['params']['pdf']['ub1 c']={'value': 0.0,'fixed':True}
conf['params']['pdf']['ub1 d']={'value': 0.0,'fixed':True}
    
conf['params']['pdf']['s1 N']={'value': 0.0,'fixed':True}
conf['params']['pdf']['s1 a']={'value': 1.35,'fixed':True}
conf['params']['pdf']['s1 b']={'value': 6.0,'fixed':True}
conf['params']['pdf']['s1 c']={'value': 0.0,'fixed':True}
conf['params']['pdf']['s1 d']={'value': 0.0,'fixed':True}
    
conf['params']['pdf']['sb1 N']={'value': 7.46,'fixed':True}
conf['params']['pdf']['sb1 a']={'value': 3.83,'fixed':True}
conf['params']['pdf']['sb1 b']={'value': 4.61,'fixed':True}
conf['params']['pdf']['sb1 c']={'value': 0.0,'fixed':True}
conf['params']['pdf']['sb1 d']={'value': 0.0,'fixed':True}
    
conf['params']['pdf']['sea1 N']={'value': 5.71,'fixed':True}
conf['params']['pdf']['sea1 a']={'value': -1.36,'fixed':True}
conf['params']['pdf']['sea1 b']={'value': 4.74,'fixed':True}
conf['params']['pdf']['sea1 c']={'value': 0.0,'fixed':True}
conf['params']['pdf']['sea1 d']={'value': 0.0,'fixed':True}
        
conf['params']['pdf']['sea2 N']={'value': 2.09,'fixed':True}
conf['params']['pdf']['sea2 a']={'value': -1.5,'fixed':True}
conf['params']['pdf']['sea2 b']={'value': 10.0,'fixed':True}
conf['params']['pdf']['sea2 c']={'value': 0.0,'fixed':True}
conf['params']['pdf']['sea2 d']={'value': 0.0,'fixed':True}

# Polarized PDFs
conf['params']['ppdf']={}
conf['params']['ppdf']['g1 N']={'value': 1.547,'fixed':True}
conf['params']['ppdf']['g1 a']={'value': -0.714,'fixed':True}
conf['params']['ppdf']['g1 b']={'value': 0.946,'fixed':True}
conf['params']['ppdf']['g1 d']={'value': 0.0,'fixed':True}
        
conf['params']['ppdf']['uv1 N']={'value': 0.348,'fixed':True}
conf['params']['ppdf']['uv1 a']={'value': 0.293,'fixed':True}
conf['params']['ppdf']['uv1 b']={'value': 2.958,'fixed':True}
conf['params']['ppdf']['uv1 d']={'value': 0.0,'fixed':True}
        
conf['params']['ppdf']['dv1 N']={'value': 0.152,'fixed':True}
conf['params']['ppdf']['dv1 a']={'value': 0.534,'fixed':True}
conf['params']['ppdf']['dv1 b']={'value': 2.011,'fixed':True}
conf['params']['ppdf']['dv1 d']={'value': 0.0,'fixed':True}
        
conf['params']['ppdf']['db1 N']={'value': -0.122,'fixed':True}
conf['params']['ppdf']['db1 a']={'value': 1.369,'fixed':True}
conf['params']['ppdf']['db1 b']={'value': 10.0,'fixed':True}
conf['params']['ppdf']['db1 d']={'value': 0.0,'fixed':True}
    
conf['params']['ppdf']['ub1 N']={'value': -0.122,'fixed':True}
conf['params']['ppdf']['ub1 a']={'value': 1.369,'fixed':True}
conf['params']['ppdf']['ub1 b']={'value': 10.0,'fixed':True}
conf['params']['ppdf']['ub1 d']={'value': 0.0,'fixed':True}
    
conf['params']['ppdf']['sea1 N']={'value': -0.122,'fixed':True}
conf['params']['ppdf']['sea1 a']={'value': 1.369,'fixed':True}
conf['params']['ppdf']['sea1 b']={'value': 10.0,'fixed':True}
conf['params']['ppdf']['sea1 d']={'value': 0.0,'fixed':True}

conf['gpd model']=2
# H GPD
conf['params']['Hgpd']={}
conf['params']['Hgpd']['g alpha_p']={'value': 0.15,'min':0.0,'max':2.0,'fixed':False,'forcemin':True,'forcemax':False}
conf['params']['Hgpd']['g b']={'value': 2.58,'min':0.0,'max':3.0,'fixed':False,'forcemin':True,'forcemax':False}
conf['params']['Hgpd']['g n']={'value':2.0,'min':0.0,'max':3.0,'fixed':False,'forcemin':True,'forcemax':False}
        
conf['params']['Hgpd']['up alpha_p']={'value': 0.9,'min':0.0,'max':2.0,'fixed':False,'forcemin':True,'forcemax':False}
conf['params']['Hgpd']['up b']={'value': 0.0,'min':0.0,'max':3.0,'fixed':False,'forcemin':True,'forcemax':False}
conf['params']['Hgpd']['up n']={'value': 1.0,'min':0.0,'max':3.0,'fixed':False,'forcemin':True,'forcemax':False}
        
conf['params']['Hgpd']['dp alpha_p']={'value': 0.9,'fixed': 'up alpha_p'}
conf['params']['Hgpd']['dp b']={'value': 0.0,'fixed': 'up b'}
conf['params']['Hgpd']['dp n']={'value': 1.0,'fixed': 'up n'}
        
conf['params']['Hgpd']['sp alpha_p']={'value': 0.15,'fixed': 'up alpha_p'}
conf['params']['Hgpd']['sp b']={'value': 2.58,'fixed': 'up b'}
conf['params']['Hgpd']['sp n']={'value': 2.0,'fixed': 'up n'}

# E GPD
conf['params']['Egpd']={}
conf['params']['Egpd']['g N']={'value': 0.388, 'min':-2.0,'max':2.0,'fixed':True}
conf['params']['Egpd']['g alpha']={'value':0.122,'min':-1.0,'max':0.999,'fixed':False, 'forcemin':False, 'forcemax':True}
conf['params']['Egpd']['g beta']={'value':3.21, 'min':0.0,'max':10.0,'fixed':False, 'forcemin':True, 'forcemax':False}
conf['params']['Egpd']['g alpha_p']={'value':0.15,'min':0.0,'max':2.0,'fixed':False,'forcemin':True,'forcemax':False}
conf['params']['Egpd']['g b']={'value':2.58,'min':0.0,'max':3.0,'fixed':False,'forcemin':True,'forcemax':False}
conf['params']['Egpd']['g n']={'value':2.0,'min':0.0,'max':3.0,'fixed':False,'forcemin':True,'forcemax':False}
        
conf['params']['Egpd']['up N']={'value':0.346, 'min':-2.0, 'max':2.0, 'fixed':False, 'forcemin':False, 'forcemax':False}
conf['params']['Egpd']['up alpha']={'value':-0.122,'min':-1.0,'max':0.999,'fixed':False, 'forcemin':False, 'forcemax':True}
conf['params']['Egpd']['up beta']={'value':3.21, 'min':0.0,'max':10.0,'fixed':False, 'forcemin':True, 'forcemax':False}
conf['params']['Egpd']['up alpha_p']={'value':0.9,'min':0.0,'max':2.0,'fixed':False,'forcemin':True,'forcemax':False}
conf['params']['Egpd']['up b']={'value':0.0,'min':0.0,'max':3.0,'fixed':False,'forcemin':True,'forcemax':False}
conf['params']['Egpd']['up n']={'value':1.0,'min':0.0,'max':3.0,'fixed':False,'forcemin':True,'forcemax':False}
        
conf['params']['Egpd']['dp N']={'value':0.152, 'fixed': 'up N'}
conf['params']['Egpd']['dp alpha']={'value':-0.24,'fixed': 'up alpha'}
conf['params']['Egpd']['dp beta']={'value':3.84,'fixed': 'up beta'}
conf['params']['Egpd']['dp alpha_p']={'value':0.9,'fixed': 'up alpha_p'}
conf['params']['Egpd']['dp b']={'value':0.0,'fixed': 'up b'}
conf['params']['Egpd']['dp n']={'value':1.0,'fixed': 'up n'}
        
conf['params']['Egpd']['sp N']={'value':0.0, 'fixed': 'up N'}
conf['params']['Egpd']['sp alpha']={'value':1.35,'fixed': 'up alpha'}
conf['params']['Egpd']['sp beta']={'value':6.0,'fixed': 'up beta'}
conf['params']['Egpd']['sp alpha_p']={'value':0.15,'fixed': 'up alpha_p'}
conf['params']['Egpd']['sp b']={'value':2.58,'fixed': 'up b'}
conf['params']['Egpd']['sp n']={'value':2.0,'fixed': 'up n'}

# Ht GPD
conf['params']['Htgpd']={}
conf['params']['Htgpd']['g alpha_p']={'value': 0.15,'min':0.0,'max':2.0,'fixed':False,'forcemin':True,'forcemax':False}
conf['params']['Htgpd']['g b']={'value': 2.58,'min':0.0,'max':3.0,'fixed':False,'forcemin':True,'forcemax':False}
conf['params']['Htgpd']['g n']={'value': 2.0,'min':0.0,'max':3.0,'fixed':False,'forcemin':True,'forcemax':False}
    
conf['params']['Htgpd']['up alpha_p']={'value': 0.9,'min':0.0,'max':2.0,'fixed':False,'forcemin':True,'forcemax':False}
conf['params']['Htgpd']['up b']={'value': 0.0,'min':0.0,'max':3.0,'fixed':False,'forcemin':True,'forcemax':False}
conf['params']['Htgpd']['up n']={'value': 1.0,'min':0.0,'max':3.0,'fixed':False,'forcemin':True,'forcemax':False}
    
conf['params']['Htgpd']['dp alpha_p']={'value': 0.9,'fixed': 'up alpha_p'}
conf['params']['Htgpd']['dp b']={'value': 0.0,'fixed': 'up b'}
conf['params']['Htgpd']['dp n']={'value': 1.0,'fixed': 'up n'}
    
conf['params']['Htgpd']['sp alpha_p']={'value': 0.15,'fixed': 'up alpha_p'}
conf['params']['Htgpd']['sp b']={'value': 2.58,'fixed': 'up b'}
conf['params']['Htgpd']['sp n']={'value': 2.0,'fixed': 'up n' }

# Et GP
conf['params']['Etgpd']={}
conf['params']['Etgpd']['g N']={'value': 0.388, 'min':-2.0, 'max':2.0, 'fixed':False, 'forcemin':False, 'forcemax':False}
conf['params']['Etgpd']['g alpha']={'value': -0.623,'min':-1.0,'max':0.999,'fixed':False, 'forcemin':False, 'forcemax':True}
conf['params']['Etgpd']['g beta']={'value': 9.26, 'min':0.0,'max':10.0,'fixed':False, 'forcemin':True, 'forcemax':False}
conf['params']['Etgpd']['g alpha_p']={'value': 0.15,'min':0.0,'max':2.0,'fixed':False,'forcemin':True,'forcemax':False}
conf['params']['Etgpd']['g b']={'value': 2.58,'min':0.0,'max':3.0,'fixed':False,'forcemin':True,'forcemax':False}
conf['params']['Etgpd']['g n']={'value': 2.0,'min':0.0,'max':3.0,'fixed':False,'forcemin':True,'forcemax':False}
        
conf['params']['Etgpd']['up N']={'value':0.346, 'min':-2.0, 'max':2.0, 'fixed':False, 'forcemin':False, 'forcemax':False}
conf['params']['Etgpd']['up alpha']={'value': -0.122,'min':-1.0,'max':0.999,'fixed':False, 'forcemin':False, 'forcemax':True}
conf['params']['Etgpd']['up beta']={'value': 3.21, 'min':0.0,'max':10.0,'fixed':False, 'forcemin':True, 'forcemax':False}
conf['params']['Etgpd']['up alpha_p']={'value': 0.9,'min':0.0,'max':2.0,'fixed':False,'forcemin':True,'forcemax':False}
conf['params']['Etgpd']['up b']={'value': 0.0,'min':0.0,'max':3.0,'fixed':False,'forcemin':True,'forcemax':False}
conf['params']['Etgpd']['up n']={'value': 1.0,'min':0.0,'max':3.0,'fixed':False,'forcemin':True,'forcemax':False}
    
conf['params']['Etgpd']['dp N']={'value': 0.152,'fixed': 'up N'}
conf['params']['Etgpd']['dp alpha']={'value': -0.24,'fixed': 'up alpha'}
conf['params']['Etgpd']['dp beta']={'value': 3.84,'fixed': 'up beta'}
conf['params']['Etgpd']['dp alpha_p']={'value': 0.9,'fixed': 'up alpha_p'}
conf['params']['Etgpd']['dp b']={'value': 0.0,'fixed': 'up b'}
conf['params']['Etgpd']['dp n']={'value': 1.0,'fixed': 'up n'}
    
conf['params']['Etgpd']['sp N']={'value': 0.0,'fixed': 'up N'}
conf['params']['Etgpd']['sp alpha']={'value': 1.35,'fixed': 'up alpha'}
conf['params']['Etgpd']['sp beta']={'value': 6.0,'fixed': 'up beta'}
conf['params']['Etgpd']['sp alpha_p']={'value': 0.15,'fixed': 'up alpha_p'}
conf['params']['Etgpd']['sp b']={'value': 2.58,'fixed': 'up b'}
conf['params']['Etgpd']['sp n']={'value': 2.0,'fixed': 'up n'}