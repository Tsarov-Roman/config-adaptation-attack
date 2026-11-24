import numpy as np, json
from sim2 import *
SYN=lam_profile(); rng=np.random.default_rng(77)
res={}
for T in [10,60]:
  for d in [0,60]:
    for P in [60,120,240,360,600,1200,1800,3600,7200]:
      if P/2<2*T+2*d: continue
      c=dict(eff=0,noint=0,cond=0,asm=0,ok=0,lo_viol=0,hi_viol=0,cens=0)
      for i in range(30):
        phi=int(rng.integers(0,P)); pinfo=np.zeros((2000,6)); pinfo[:,5]=-1
        r=run2(SYN,'R3',T=T,d=d,amode=1,a=250,P=P,phi=phi,seed=50000+i,pinfo=pinfo)
        pi=pinfo[:r[6]]
        for th,v,cen,itm,eff,td in pi:
            if eff!=1: continue
            c['eff']+=1
            if cen==1: c['cens']+=1; continue
            if itm!=1: continue
            c['noint']+=1
            lastfull = T-1 if th==0 else 2*T-1-th
            cond = 0<=td<=lastfull
            if cond: c['cond']+=1
            lo=T-th+d; hi=2*T-th+d
            if cond:
                c['asm']+=1
                if lo<=v<=hi: c['ok']+=1
                else:
                    if v<lo: c['lo_viol']+=1
                    if v>hi: c['hi_viol']+=1
      res[f'T{T}|d{d}|P{P}']=c
json.dump(res,open('verify_p2.json','w'))
tot={k:sum(r[k] for r in res.values()) for k in ['eff','cens','noint','cond','asm','ok','lo_viol','hi_viol']}
print(tot)
fr=[r['cond']/r['noint'] for r in res.values()]; print('cond share among noint: %.4f..%.4f'%(min(fr),max(fr)))
print('ok==asm everywhere:', all(r['ok']==r['asm'] for r in res.values()))
