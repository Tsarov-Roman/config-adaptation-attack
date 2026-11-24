import numpy as np, json, time
from scipy import stats
from sim2 import *
from model import opt_index, CONF
from trace import days
t0=time.time()
SYN=lam_profile(); TR,TRN=days()
OUT={}
def ci(v):
    v=np.asarray(v,float); return [float(v.mean()), float(stats.t.ppf(0.975,len(v)-1)*v.std(ddof=1)/np.sqrt(len(v))) if v.std()>0 else 0.0]
# ---- E7: per-pulse check of Prop 2 (R3)
rng=np.random.default_rng(77)
E7={}
for T in [10,60]:
  for d in [0,60]:
    for P in [60,120,240,360,600,1200,1800,3600,7200]:
      if P/2<2*T+2*d: continue
      tinit=[]; q=[]; tot=[]; nin=0; nok=0; nasm=0; neff=0; nall=0
      for i in range(30):
        phi=int(rng.integers(0,P)); pinfo=np.zeros((2000,5))
        r=run2(SYN,'R3',T=T,d=d,amode=1,a=250,P=P,phi=phi,seed=50000+i,pinfo=pinfo)
        npl=r[6]; pi=pinfo[:npl]
        eff=pi[:,4]==1; neff+=eff.sum(); nall+=npl
        v=pi[eff,1]; th=pi[eff,0]; cens=pi[eff,2]==1; interm=pi[eff,3]
        tinit.append(v.sum()/DAY); q.append(eff.mean()); tot.append(r[1])
        asm=(~cens)&(interm==1); nasm+=asm.sum()
        ok=(v>=T-th+d)&(v<=2*T-th+d); nok+=(ok&asm).sum()
      E7[f'T{T}|d{d}|P{P}']={'tinit':ci(np.array(tinit)*100),'tot':ci(np.array(tot)*100),'q':float(np.mean(q)),
         'share_asm':float(nasm/max(neff,1)),'share_ok_given_asm':float(nok/max(nasm,1)),'lo':100*(T/2+d)/P,'hi':100*(1.5*T+d)/P}
OUT['E7']=E7; print('E7',time.time()-t0)
# ---- E8: noisy boundary-tracking attacker
E8={}
for kind,prof in [('syn',[SYN]*100),('tr',TR)]:
  for sr in [0.0,0.005,0.01,0.02,0.05]:
    for name,kw in [('R4',{}),('PR',{}),('HPA',{}),('R5',{'ares':250}),('PR5',{'ares':250})]:
      R=[]
      for i,l in enumerate(prof):
        cb=base_cost(l); r=run2(l,name,amode=3,a=250,sig_rel=sr,obs_dt=10,seed=60000+i,**kw)
        R.append((100*r[1],100*r[3]/l.mean(),100*(r[2]/cb-1),100*r[4],100*r[5]))
      R=np.array(R)
      E8[f'{kind}|{sr}|{name}']={'viol':ci(R[:,0]),'share':ci(R[:,1]),'cost':ci(R[:,2]),'on':ci(R[:,3]),'onv':ci(R[:,4]),
        'prec':float(R[:,4].sum()/max(R[:,3].sum(),1e-9)),'raw_viol':R[:,0].tolist()}
OUT['E8']=E8; print('E8',time.time()-t0)
# ---- E9: lower margin rho_minus
xs=np.linspace(1,6000,600001); ii=opt_index(xs)
brk=xs[1:][ii[1:]!=ii[:-1]]  # first lam of new config
def rho_minus(lam):
    k=np.searchsorted(brk,lam,side='right')-1
    L=np.where(k>=0,brk[np.clip(k,0,None)],0.0)
    return (lam-L)/lam
def rho_plus(lam):
    j=opt_index(lam); return (CONF[j,3]-lam)/lam
E9={}
for nm,prof in [('syn',[SYN]),('tr',TR)]:
    rm=np.concatenate([rho_minus(l[::10]) for l in prof]); rp=np.concatenate([rho_plus(l[::10]) for l in prof])
    E9[nm]={'minus':{'median':float(np.median(rm)),'p10':float(np.quantile(rm,0.1)),'p5':float(np.quantile(rm,0.05)),'lt1':float(np.mean(rm<0.01)),'lt2':float(np.mean(rm<0.02))},
            'plus':{'median':float(np.median(rp)),'p10':float(np.quantile(rp,0.1)),'p5':float(np.quantile(rp,0.05))}}
OUT['E9']=E9
np.save('rho_syn.npy',rho_minus(SYN)); np.save('rho_tr.npy',np.concatenate([rho_minus(l[::10]) for l in TR]))
# ---- E10: raw trace data for block bootstrap (same seeds as experiments.py -> reproduce)
json.dump(OUT,open('results_r4.json','w')); print('done',time.time()-t0)
