"""Additional runs used in the paper: PR+reserve (results_pr5.json) and the reserve cost table (results_t4.json).
The code below reproduces the two inline scripts that produced these files."""
import numpy as np, json
from scipy import stats
from sim import *
from model import opt_index
from trace import days
SYN=lam_profile(); TR,_=days()
def ci(v):
    v=np.asarray(v,float); return [float(v.mean()), float(stats.t.ppf(0.975,len(v)-1)*v.std(ddof=1)/np.sqrt(len(v))) if v.std()>0 else 0.0]
def pr5():
    rng=np.random.default_rng(7); out={}
    def rs(prof,rule,**kw):
        R=[]
        for i,l in enumerate(prof):
            cb=base_cost(l); n,v,c,s=run(l,rule,phi=int(rng.integers(0,600)),seed=30000+i,**kw)
            R.append((n,100*v,100*(c/cb-1),100*s/l.mean()))
        R=np.array(R); return {k:ci(R[:,j]) for j,k in enumerate(['rec','viol','cost','share'])}
    for kind,prof in [('syn',[SYN]*100),('tr',TR)]:
        for d in [0,60]:
            out[f'{kind}|d{d}|none']=rs(prof,'PR5',ares=250,d=d)
            out[f'{kind}|d{d}|P600']=rs(prof,'PR5',ares=250,d=d,amode=1,a=250,P=600)
            out[f'{kind}|d{d}|adapt']=rs(prof,'PR5',ares=250,d=d,amode=2,a=250)
    json.dump(out,open('results_pr5.json','w'))
def t4():
    rng=np.random.default_rng(11); out={}
    for kind,prof in [('syn',[SYN]*100),('tr',TR)]:
        for a in [100,250,500]:
            rn=[];rp=[];idl=[]
            for i,l in enumerate(prof):
                cb=base_cost(l); phi=int(rng.integers(0,600))
                _,v,c,_=run(l,'R5',ares=a,seed=40000+i); rn.append(100*(c/cb-1))
                _,v,c,_=run(l,'R5',ares=a,seed=40000+i,amode=1,a=a,P=600,phi=phi); rp.append(100*(c/cb-1))
                t=np.arange(len(l)); A=a*(((t+phi)%600)<300); idl.append(100*(CC[opt_index(l+A)].mean()/cb-1))
            out[f'{kind}|a{a}']={'R5none':ci(rn),'R5per':ci(rp),'ideal':ci(idl),'share':a/2/1450*100}
    json.dump(out,open('results_t4.json','w'))
if __name__=='__main__':
    pr5(); t4()
