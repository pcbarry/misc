
def _hook10(params,order):
    sample=[]
    for i in range(len(order)):
        if order[i][0]!=1: continue
        #if order[i][1]!='pdf': continue
        #if 'g1'  in order[i][2]: continue 
        #if 'uv1' in order[i][2]: continue 
        #if 'dv1' in order[i][2]: continue 
        #if 'db1' in order[i][2]: continue 
        #if 'ub1' in order[i][2]: continue 
        if 's1' in order[i][2]: continue 
        #if 'sb1' in order[i][2]: continue 
        sample.append(params[i])
    return sample

def _hook11(params,order):
    sample=[]
    for i in range(len(order)):
        if order[i][0]!=1: continue
        #if order[i][1]!='pdf': continue
        if 'g1'  in order[i][2]: continue 
        #if 'uv1' in order[i][2]: continue 
        #if 'dv1' in order[i][2]: continue 
        #if 'db1' in order[i][2]: continue 
        #if 'ub1' in order[i][2]: continue 
        #if 's1' in order[i][2]: continue 
        #if 'sb1' in order[i][2]: continue 
        sample.append(params[i])
    return sample

def _hook13(params,order):
    sample=[]
    for i in range(len(order)):
        if order[i][0]!=1: continue
        if order[i][1]!='pdf-pion': continue
        if 'ubv1' in order[i][2]: continue 
        if 'u1' in order[i][2]: continue 
        sample.append(params[i])
    return sample

nc={}
hooks={}

nc[10]=2
hooks[10]=_hook10

nc[11]=1
hooks[11]=None#_hook11

nc[12]=1
hooks[12]=None#_hook11

nc[13]=1
hooks[13]=None#_hook13

nc[14]=1
hooks[14]=None#_hook11

nc[15]=1
hooks[15]=None#_hook11

nc[16]=1
hooks[16]=None#_hook11

nc[17]=1
hooks[17]=None#_hook11
