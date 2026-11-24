import json,numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from sim import *
from model import opt_index
from trace import days
plt.rcParams.update({'font.family':'Times New Roman','font.size':9,'axes.linewidth':0.6,'lines.linewidth':0.9,
  'xtick.direction':'in','ytick.direction':'in','savefig.dpi':300,'mathtext.fontset':'stix'})
O=json.load(open('results.json')); OUT='../fig/'
lam=lam_profile()
# Fig1 scheme
fig,ax=plt.subplots(figsize=(6.3,2.9)); ax.axis('off'); ax.set_xlim(0,100); ax.set_ylim(0,47)
def box(x,y,w,h,t):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.4',fc='white',ec='black',lw=0.7))
    ax.text(x+w/2,y+h/2,t,ha='center',va='center',fontsize=8.5)
def arr(x1,y1,x2,y2,ls='-',rad=0.0):
    ax.annotate('',xy=(x2,y2),xytext=(x1,y1),arrowprops=dict(arrowstyle='-|>',lw=0.7,ls=ls,color='black',connectionstyle=f'arc3,rad={rad}'))
box(1,22,19,9,'Легітимне\nнавантаження λ(t)')
box(31,22,18,9,'Платформа,\nконфігурація x(t)')
box(60,22,17,9,'Моніторинг:\nоцінка $\\hat\\lambda$, σ')
box(60,3,17,9,'Правило\nоновлення')
box(31,3,18,9,'Реконфігурація\nтривалістю d')
box(84,22,15,9,'Нападник:\nімпульси A(t)')
box(84,3,15,9,'Нападник:\nзміщення b')
arr(20,26.5,31,26.5); arr(49,26.5,60,26.5); arr(68.5,22,68.5,12); arr(60,7.5,49,7.5); arr(40,12,40,22)
arr(89,31.6,44,31.6,'--',0.32); ax.text(66.5,35.0,'паразитні запити',fontsize=7.5,ha='center')
arr(37,31.6,95,31.6,':',-0.32); ax.text(66.5,44.3,'спостереження x(t)',fontsize=7.5,ha='center')
arr(84,9,77,22,'--')
fig.savefig(OUT+'fig1.png',bbox_inches='tight'); plt.close(fig)
# Fig2 time series
t0,t1=10*3600,12*3600
fig,axs=plt.subplots(3,1,figsize=(6.3,4.8),sharex=True)
cases=[('R3',dict(amode=1,a=250,P=600,phi=100),'а) R3, періодичні імпульси, P = 600 с'),
       ('R4',dict(amode=2,a=250),'б) R4, нападник, що відстежує межу'),
       ('R5',dict(amode=1,a=250,P=600,phi=100,ares=250),'в) R5 із запасом a$_{рез}$ = a, періодичні імпульси')]
for ax,(r,kw,tt) in zip(axs,cases):
    rec=np.zeros((2,DAY)); run(lam,r,seed=5,recarr=rec,**kw)
    th=np.arange(DAY)/3600
    ax.plot(th[t0:t1],rec[0,t0:t1],color='0.45',lw=0.6,label='λ + A')
    ax.step(th[t0:t1],rec[1,t0:t1],where='post',color='black',lw=0.9,label='λ$_{max}$(x)')
    v=rec[0]>rec[1]
    ax.fill_between(th[t0:t1],1500,2700,where=v[t0:t1],color='0.85',step='post',lw=0)
    ax.set_ylim(1500,2700); ax.set_ylabel('запит/с'); ax.set_title(tt,fontsize=9,loc='left')
axs[0].legend(loc='upper left',fontsize=8,frameon=False,ncol=2); axs[-1].set_xlabel('час доби, год')
fig.tight_layout(); fig.savefig(OUT+'fig2.png'); plt.close(fig)
# Fig3 period sweep
fig,axs=plt.subplots(1,2,figsize=(6.3,2.6),sharey=True)
Ps=np.array([60,120,240,360,600,1200,1800,3600,7200])
for ax,T in zip(axs,[10,60]):
    for d,mk in [(0,'o'),(60,'s')]:
        v=[O[f'E3|T{T}|d{d}|P{P}|R3']['viol'][0] for P in Ps]
        v1=[O[f'E3|T{T}|d{d}|P{P}|R1']['viol'][0] for P in Ps]
        ax.loglog(Ps,v,mk+'-',color='black',ms=3.5,label=f'R3, d = {d} с')
        ax.loglog(Ps,v1,mk+':',color='0.5',ms=3,mfc='white',label=f'R1, d = {d} с')
        ok=Ps/2>=2*T+2*d; pp=Ps[ok]
        ax.fill_between(pp,100*(T/2+d)/pp,100*(1.5*T+2*d)/pp,color='0.88' if d==0 else '0.78',alpha=0.7,lw=0)
    ax.set_title(f'T = {T} с',fontsize=9); ax.set_xlabel('період P, с'); ax.set_ylim(0.05,60)
axs[0].set_ylabel('час порушення, %'); axs[1].legend(fontsize=7,frameon=False,loc='lower left')
fig.tight_layout(); fig.savefig(OUT+'fig3.png'); plt.close(fig)
# Fig4 cost-violation plane
fig,ax=plt.subplots(figsize=(6.3,3.0))
names=['R1','R2','R3','R4','R5','PR','HPA']
for nm in names:
    p=O[f'E2|syn|d0|P600|{nm}']; q=O.get(f'E4|syn|d0|a250|{nm}')
    off={'R5':(-14,-9),'PR':(4,4)}.get(nm,(3,2))
    ax.plot(p['cost'][0],p['viol'][0],'o',color='black',ms=4); ax.annotate(nm,(p['cost'][0],p['viol'][0]),xytext=off,textcoords='offset points',fontsize=7.5)
    if q:
        ax.plot(q['cost'][0],q['viol'][0],'^',color='black',mfc='white',ms=4.5); ax.annotate(nm,(q['cost'][0],q['viol'][0]),xytext=((-14,-9) if nm=='R5' else (3,-8)),textcoords='offset points',fontsize=7.5,color='0.35')
pr5=json.load(open('results_pr5.json'))
for k,mk in [('syn|d0|P600','o'),('syn|d0|adapt','^')]:
    ax.plot(pr5[k]['cost'][0],pr5[k]['viol'][0],mk,color='black',mfc='black' if mk=='o' else 'white',ms=4.5)
    ax.annotate('PR+рез',(pr5[k]['cost'][0],pr5[k]['viol'][0]),xytext=(3,3),textcoords='offset points',fontsize=7.5)
lb=O['E5|syn|a250'][0]; ideal=O['E2|syn|P600|ideal'][0]
ax.axvline(lb,ls='--',color='black',lw=0.7); ax.text(lb+0.4,45,'нижня межа\n(твердження 3)',fontsize=7.5)
ax.axvline(ideal,ls=':',color='0.4',lw=0.7); ax.text(ideal+0.4,52,'ідеальна реакція',fontsize=7.5,color='0.3')
ax.set_xlabel('приріст вартості, %'); ax.set_ylabel('час порушення, %'); ax.set_ylim(-2,56)
ax.plot([],[],'o',color='black',ms=4,label='періодичні імпульси, P = 600 с'); ax.plot([],[],'^',color='black',mfc='white',ms=4.5,label='нападник, що відстежує межу')
ax.legend(fontsize=7.5,frameon=False,loc='upper right')
fig.tight_layout(); fig.savefig(OUT+'fig4.png'); plt.close(fig)
# Fig5 rho + bias
fig,axs=plt.subplots(1,2,figsize=(6.3,2.6))
TR,_=days()
r_s=(LM[opt_index(lam)]-lam)/lam*100
r_t=np.concatenate([(LM[opt_index(l[::10])]-l[::10])/l[::10]*100 for l in TR])
for r,ls,lb_ in [(r_s,'-','синтетичний профіль'),(r_t,'--','трасa WC98, 43 доби')]:
    xs=np.sort(r); axs[0].plot(xs,np.arange(len(xs))/len(xs),ls,color='black',label=lb_.replace('трасa','траса'))
axs[0].set_xlim(0,16); axs[0].set_xlabel('запас ρ, %'); axs[0].set_ylabel('частка доби'); axs[0].legend(fontsize=7.5,frameon=False,loc='lower right')
axs[0].set_title('а)',fontsize=9,loc='left')
bs=[-0.10,-0.05,-0.03,-0.02,-0.01,0.0]
for kind,ls in [('syn','-'),('tr','--')]:
    ex=[O[f'E6|{kind}|b{b}']['viol_exact'][0] for b in bs]; r3=[O[f'E6|{kind}|b{b}']['R3']['viol'][0] for b in bs]
    axs[1].plot([-100*b for b in bs],ex,'o'+ls,color='black',ms=3,label=('точна оцінка' if kind=='syn' else None))
    axs[1].plot([-100*b for b in bs],r3,'s'+ls,color='0.5',ms=3,mfc='white',label=('R3' if kind=='syn' else None))
axs[1].set_xlabel('заниження оцінки |b|, %'); axs[1].set_ylabel('час порушення, %'); axs[1].plot([],[],'-',color='black',label='синтетичний профіль'); axs[1].plot([],[],'--',color='black',label='траса WC98'); axs[1].legend(fontsize=7.5,frameon=False)
axs[1].set_title('б)',fontsize=9,loc='left')
fig.tight_layout(); fig.savefig(OUT+'fig5.png'); plt.close(fig)
print('ok')
