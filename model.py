import numpy as np
from scipy.stats import binom
MUS=[100,150,200,300,400,600]; A_COST=40; B_COST=0.005; C_CH=60
TAU=0.02; TMAX=0.05; MTBF=2000; MTTR=8; AS=MTBF/(MTBF+MTTR); PREQ=0.9999; NMAX=60
def cost(N,mu): return N*(A_COST+B_COST*mu**1.8)+C_CH*np.ceil(N/4)
def lam_max(N,mu):
    lrt=N*AS*(mu-1/(TMAX-TAU))
    # largest m with P(X>=m)>=PREQ
    mm=0
    for m in range(N,-1,-1):
        if binom.sf(m-1,N,AS)>=PREQ: mm=m; break
    lav=0.95*mu*mm
    return max(0.0,min(lrt,lav))
CONF=[(N,mu,cost(N,mu),lam_max(N,mu)) for mu in MUS for N in range(1,NMAX+1)]
CONF=np.array(CONF)
# tie-break fixed order: by cost then by index
order=np.lexsort((np.arange(len(CONF)),CONF[:,2]))
CONF=CONF[order]
# envelope: for lam, cheapest config with lammax>=lam
idx=np.argsort(-CONF[:,3],kind='stable')
lm_sorted=CONF[idx,3]
best=[];cur=None
for i in idx:
    if cur is None or CONF[i,2]<CONF[cur,2]: cur=i
    best.append(cur)
best=np.array(best)
def opt_index(lam):
    lam=np.atleast_1d(lam)
    # number of configs with lammax>=lam -> position
    k=np.searchsorted(-lm_sorted,-lam,side='right')-1
    out=np.where(k>=0,best[np.clip(k,0,None)],-1)
    return out
def seq(lo,hi):
    xs=np.linspace(lo,hi,200001); ii=opt_index(xs)
    ch=np.r_[True,ii[1:]!=ii[:-1]]
    return [(int(CONF[i,0]),int(CONF[i,1]),round(CONF[i,3],1)) for i in ii[ch]]
if __name__=='__main__':
    s=seq(600,3486); print(len(s)); print(s)
