import json,numpy as np
O=json.load(open('results_full.json'))
rng=np.random.default_rng(2027)
def mbb_ci(x,B=10000,L=7):
    x=np.asarray(x,float); n=len(x); nb=int(np.ceil(n/L)); st=rng.integers(0,n-L+1,(B,nb))
    idx=(st[:,:,None]+np.arange(L)).reshape(B,-1)[:,:n]; m=x[idx].mean(1); return np.percentile(m,[2.5,97.5])
def ac1(x): x=np.asarray(x); return float(np.corrcoef(x[:-1],x[1:])[0,1])
S={}
H1=np.array(O['E2|tr|d0|PNone|R3']['raw_viol'])-np.array(O['E2|tr|d0|PNone|PR']['raw_viol'])
H4=np.array(O['E4|tr|d0|a250|R5']['raw_viol'])-np.array(O['E2|tr|d0|P600|R5']['raw_viol'])
for nm,d in [('H1',H1),('H4',H4)]:
    S[nm]={'mean':float(d.mean()),'median':float(np.median(d)),'ac1':ac1(d),'ci':mbb_ci(d).tolist()}
# block half-widths for every trace table cell
HW={}
for k,v in O.items():
    if isinstance(v,dict) and '|tr|' in k:
        for f in ('raw_viol','raw_cost','raw_rec','raw_share'):
            if f in v:
                lo,hi=mbb_ci(v[f]); HW[k+'|'+f]=[float(np.mean(v[f])),float((hi-lo)/2)]
        if 'raw_exv' in v:
            lo,hi=mbb_ci(v['raw_exv']); HW[k+'|raw_exv']=[float(np.mean(v['raw_exv'])),float((hi-lo)/2)]
        if 'R3' in v:
            for f in ('raw_viol','raw_cost'):
                lo,hi=mbb_ci(v['R3'][f]); HW[k+'|R3|'+f]=[float(np.mean(v['R3'][f])),float((hi-lo)/2)]
json.dump({'S':S,'HW':HW},open('stats_r4.json','w'))
print(S)
for k in ['E2|tr|d0|PNone|R3|raw_viol','E2|tr|d0|PNone|PR|raw_viol','E4|tr|d0|a250|R5|raw_viol']: print(k,HW[k])
