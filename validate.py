import numpy as np
from sim import *
lam=lam_profile(); cb=base_cost(lam)
rng=np.random.default_rng(1)
def agg(rule,**kw):
    r=[]
    for i in range(20):
        P=kw.get('P',600); phi=int(rng.integers(0,P))
        n,v,c,_=run(lam,rule,phi=phi,seed=1000+i,**kw); r.append((n,v*100,(c/cb-1)*100))
    return np.mean(r,0).round(2)
for R in ['R1','R2','R3']: print('noattack',R,agg(R))
for P in [120,600,3600]:
    for R,kw in [('R1',{}),('R2',{}),('R3',{}),('R4',{}),('R5',{'ares':250}),('R5h',{'ares':125})]:
        print(P,R,agg(R[:2],amode=1,a=250,P=P,**kw))
