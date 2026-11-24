import numpy as np, json, time
from scipy import stats
from sim import *
from model import opt_index
from trace import days
t0=time.time()
rng=np.random.default_rng(2026)
SYN=lam_profile(); TR,TRN=days()
def ci(v):
    v=np.asarray(v,float); m=v.mean()
    if len(v)<2 or v.std()==0: return [m,0.0]
    return [m, stats.t.ppf(0.975,len(v)-1)*v.std(ddof=1)/np.sqrt(len(v))]
def runset(profiles,rule,phases,seeds,**kw):
    R=[]
    for lam,phi,sd in zip(profiles,phases,seeds):
        cb=base_cost(lam)
        n,v,c,s=run(lam,rule,phi=int(phi),seed=int(sd),**kw)
        R.append((n,100*v,100*(c/cb-1),100*s/lam.mean()))
    return np.array(R)
def summ(R): return {'rec':ci(R[:,0]),'viol':ci(R[:,1]),'cost':ci(R[:,2]),'share':ci(R[:,3]),'raw_viol':R[:,1].tolist(),'raw_cost':R[:,2].tolist()}
NS=100
def setup(kind,P=600):
    if kind=='syn':
        prof=[SYN]*NS; ph=rng.integers(0,10**6,NS)%P; sd=np.arange(NS)+10000
    else:
        prof=TR; n=len(TR); ph=rng.integers(0,10**6,n)%P; sd=np.arange(n)+20000
    return prof,ph,sd
RULES=[('R1',{}),('R2',{}),('R3',{}),('R4',{}),('R5',{'ares':'a'}),('R5h',{'ares':'a2'}),('PR',{}),('HPA',{})]
def kwr(name,kw,a):
    k=dict(kw)
    if 'ares' in k: k['ares']=a if k['ares']=='a' else a/2
    if name=='R5h': name='R5'
    return name,k
OUT={}
# E0 configurations, margin
xs=np.linspace(600,3486,200001); ii=opt_index(xs)
OUT['n_conf_600_3486']=int(len(np.unique(ii)))
def rho(lam):
    j=opt_index(lam); return (LM[j]-lam)/lam
for nm,prof in [('syn',[SYN]),('tr',TR)]:
    r=np.concatenate([rho(l[::10]) for l in prof])
    OUT['rho_'+nm]={'median':float(np.median(r)),'p10':float(np.quantile(r,0.1)),'p5':float(np.quantile(r,0.05)),'hist':np.histogram(r*100,bins=np.arange(0,16.5,0.5))[0].tolist()}
# E1/E2 rules x attack x delay
for kind in ['syn','tr']:
    for d in [0,60,300]:
        for P in [None,120,600,3600]:
            if kind=='tr' and P not in (None,600): continue
            prof,ph,sd=setup(kind,P or 600)
            for name,kw in RULES:
                if P is None and name=='R5h': continue
                rn,k=kwr(name,kw,250)
                if P is None:
                    R=runset(prof,rn,ph,sd,d=d,**k)
                else:
                    R=runset(prof,rn,ph,sd,d=d,amode=1,a=250,P=P,**k)
                OUT[f'E2|{kind}|d{d}|P{P}|{name}']=summ(R)
            if P:
                # ideal response
                ic=[]
                for lam,phi in zip(prof,ph):
                    t=np.arange(len(lam)); A=250*(((t+phi)%P)<P/2)
                    ic.append(100*(CC[opt_index(lam+A)].mean()/base_cost(lam)-1))
                OUT[f'E2|{kind}|P{P}|ideal']=ci(ic)
print('E2 done',time.time()-t0)
# E3 period sweep R3, R1
for T in [10,60]:
    for d in [0,60]:
        for P in [60,120,240,360,600,1200,1800,3600,7200]:
            prof,ph,sd=setup('syn',P); prof,ph,sd=prof[:30],ph[:30],sd[:30]
            for rn in ['R1','R3']:
                R=runset(prof,rn,ph,sd,T=T,d=d,amode=1,a=250,P=P)
                OUT[f'E3|T{T}|d{d}|P{P}|{rn}']={'viol':ci(R[:,1]),'rec':ci(R[:,0])}
print('E3 done',time.time()-t0)
# E4 adaptive attacker
for kind in ['syn','tr']:
    for d in [0,60]:
        for a in [100,250,500]:
            prof,ph,sd=setup(kind)
            for name,kw in [('R3',{}),('R4',{}),('R5',{'ares':'a'}),('PR',{}),('HPA',{})]:
                rn,k=kwr(name,kw,a)
                R=runset(prof,rn,ph,sd,d=d,amode=2,a=a,theta=0,**k)
                OUT[f'E4|{kind}|d{d}|a{a}|{name}']=summ(R)
                if a==250:
                    Rp=runset(prof,rn,ph,sd,d=d,amode=1,a=a,P=600,**k)
                    # efficiency: violation per share
                    ea=R[:,1]/np.maximum(R[:,3],1e-9); ep=Rp[:,1]/Rp[:,3]
                    w=stats.wilcoxon(ea,ep) if np.any(ea!=ep) else None
                    OUT[f'E4eff|{kind}|d{d}|{name}']={'adaptive':ci(ea),'periodic':ci(ep),'p':(float(w.pvalue) if w else None)}
# theta sensitivity
prof,ph,sd=setup('syn')
for th in [0,15,30,45,55]:
    for name in ['R3','R4','PR','HPA']:
        R=runset(prof[:30],name,ph[:30],sd[:30],amode=2,a=250,theta=th)
        OUT[f'E4th|{th}|{name}']={'viol':ci(R[:,1]),'share':ci(R[:,3])}
print('E4 done',time.time()-t0)
# E5 lower bound for guaranteed admissibility
for kind,prof in [('syn',[SYN]),('tr',TR)]:
    for a in [100,250,500]:
        lb=[100*(CC[opt_index(l+a)].mean()/base_cost(l)-1) for l in prof]
        OUT[f'E5|{kind}|a{a}']=ci(lb)
# E6 bias
for kind in ['syn','tr']:
    prof,ph,sd=setup(kind)
    for b in [-0.10,-0.05,-0.03,-0.02,-0.01,0.0,0.05,0.10]:
        ex=[];exv=[]
        for l in (prof[:1] if kind=='syn' else prof):
            j0=opt_index(l); j1=opt_index(l*(1+b))
            ex.append(100*np.mean(j0!=j1)); exv.append(100*np.mean(l>LM[j1]))
        R=runset(prof,'R3',ph,sd,bias=b)
        OUT[f'E6|{kind}|b{b}']={'changed':ci(ex),'viol_exact':ci(exv),'R3':summ(R)}
# paired tests key
def wil(k1,k2,f='raw_viol'):
    a=np.array(OUT[k1][f]); b=np.array(OUT[k2][f])
    return float(stats.wilcoxon(a,b).pvalue) if np.any(a!=b) else None
OUT['tests']={
 'syn_d0_P600_PR_vs_R3_viol':wil('E2|syn|d0|P600|PR','E2|syn|d0|P600|R3'),
 'syn_d0_P600_PR_vs_R5_cost':wil('E2|syn|d0|P600|PR','E2|syn|d0|P600|R5','raw_cost'),
 'syn_d0_P600_R4_vs_R3_viol':wil('E2|syn|d0|P600|R4','E2|syn|d0|P600|R3'),
 'syn_d0_P600_HPA_vs_R3_viol':wil('E2|syn|d0|P600|HPA','E2|syn|d0|P600|R3'),
 'tr_d0_P600_PR_vs_R3_viol':wil('E2|tr|d0|P600|PR','E2|tr|d0|P600|R3'),
}
for k in list(OUT):
    if isinstance(OUT[k],dict):
        OUT[k].pop('raw_viol',None) if False else None
json.dump(OUT,open('results.json','w'))
print('done',time.time()-t0)
