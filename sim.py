import numpy as np
from numba import njit
from model import CONF, lm_sorted, best, opt_index
LM=CONF[:,3].copy(); CC=CONF[:,2].copy()
LMS=lm_sorted.copy(); BEST=best.astype(np.int64)
DAY=86400

@njit(cache=True)
def opt(lam,LMS,BEST):
    # count of configs with lammax>=lam, LMS descending
    lo=0; hi=LMS.shape[0]
    while lo<hi:
        m=(lo+hi)//2
        if LMS[m]>=lam: lo=m+1
        else: hi=m
    k=lo-1
    if k<0: return -1
    return BEST[k]

@njit(cache=True)
def simulate(lam,T,rule,z,h,H,ares,d,bias,amode,a,P,phi,theta,seed,
             LM,CC,LMS,BEST,alpha,beta,mhpa,S,recarr):
    np.random.seed(seed)
    n=lam.shape[0]
    cur=opt(lam[0],LMS,BEST); pend=-1; tsw=-1
    nrec=0; viol=0; csum=0.0; asum=0.0
    acc=0.0; lastup=-10**9
    c1=-1; c2=-1
    Lv=lam[0]; Bv=0.0; sres=0.0; fprev=-1.0
    nb=max(1,S//T); ring=np.zeros(nb); ri=0; rfill=0
    on=False
    for t in range(n):
        if pend>=0 and t>=tsw:
            cur=pend; pend=-1; nrec+=1
        # attack
        A=0.0
        if amode==1:
            if ((t+phi)%P)<P/2: A=a
        elif amode==2:
            if on:
                if lam[t]+a<=LM[cur]: on=False
            else:
                if (t%T)==theta and lam[t]+a>LM[cur] and pend<0: on=True
            if on: A=a
        L=lam[t]+A; asum+=A
        if recarr.shape[1]==n:
            recarr[0,t]=L; recarr[1,t]=LM[cur]
        if L>LM[cur]: viol+=1
        c=CC[cur]
        if pend>=0: c+=CC[pend]
        csum+=c
        acc+=L
        if (t+1)%T==0:
            cnt=np.random.poisson(acc); acc=0.0
            lh=cnt/T*(1.0+bias)
            sig=np.sqrt(max(lh,1.0)/T)
            x=pend if pend>=0 else cur
            tg=-1
            if rule==0:
                tg=opt(lh,LMS,BEST)
            elif rule==1:
                cc=opt(lh,LMS,BEST)
                if cc==c1 and cc==c2: tg=cc
                c2=c1; c1=cc
            elif rule==2 or rule==3:
                if rule==2:
                    u=lh+z*sig+ares
                else:
                    if fprev>0: sres=0.9*sres+0.1*abs(lh-fprev)
                    Ln=alpha*lh+(1-alpha)*(Lv+Bv)
                    Bv=beta*(Ln-Lv)+(1-beta)*Bv; Lv=Ln
                    k=1.0+d/T
                    f=Lv+Bv*k; fprev=Lv+Bv
                    u=max(f,lh)+z*max(sig,1.25*sres)+ares
                if u>LM[x]:
                    tg=opt(u,LMS,BEST); lastup=t
                elif t-lastup>=H:
                    cc=opt((1+h)*u,LMS,BEST)
                    if CC[cc]<CC[x]: tg=cc
            elif rule==4:
                ring[ri]=lh*(1+mhpa); ri=(ri+1)%nb
                if rfill<nb: rfill+=1
                u=0.0
                for j in range(rfill):
                    if ring[j]>u: u=ring[j]
                tg=opt(u,LMS,BEST)
            if pend<0 and tg>=0 and tg!=cur:
                pend=tg; tsw=t+1+d
    return nrec, viol/n, csum/n, asum/n

def lam_profile():
    t=np.arange(DAY)+0.5
    return 1450-750*np.cos(2*np.pi*(t-4*3600)/DAY)

RULES={'R1':(0,0,0.0),'R2':(1,0,0.0),'R3':(2,0,0.0),'R4':(2,600,0.0)}
def run(lam,rule='R3',T=60,z=2.0,h=0.02,H=0,ares=0.0,d=0,bias=0.0,amode=0,a=0.0,P=600,phi=0,theta=0,seed=0,
        alpha=0.5,beta=0.2,mhpa=0.1,S=300,recarr=None):
    if recarr is None: recarr=np.zeros((2,1))
    code={'R1':0,'R2':1,'R3':2,'R4':2,'R5':2,'PR':3,'PR5':3,'HPA':4}[rule]
    if rule=='R4' and H==0: H=600
    return simulate(lam,T,code,z,h,H,ares,d,bias,amode,a,P,phi,theta,seed,LM,CC,LMS,BEST,alpha,beta,mhpa,S,recarr)

def base_cost(lam):
    return CC[opt_index(lam)].mean()
