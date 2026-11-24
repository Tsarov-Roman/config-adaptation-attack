import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sim import *
plt.rcParams.update({'font.family':'Times New Roman','font.size':9,'axes.linewidth':0.6,'lines.linewidth':0.9,
  'xtick.direction':'in','ytick.direction':'in','savefig.dpi':300,'mathtext.fontset':'stix'})
lam=lam_profile(); t0,t1=10*3600,12*3600
fig,axs=plt.subplots(3,1,figsize=(6.3,4.8),sharex=True)
cases=[('R3',dict(amode=1,a=250,P=600,phi=100),'а) R3, періодичні імпульси, P = 600 с'),
       ('R5',dict(amode=1,a=250,P=600,phi=100,ares=250),'б) R5 з резервом a$_{рез}$ = a, періодичні імпульси'),
       ('R4',dict(amode=2,a=250),'в) R4, нападник, що відстежує межу')]
for ax,(r,kw,tt) in zip(axs,cases):
    rec=np.zeros((2,DAY)); run(lam,r,seed=5,recarr=rec,**kw)
    th=np.arange(DAY)/3600
    ax.plot(th[t0:t1],rec[0,t0:t1],color='0.45',lw=0.6,label='λ + A')
    ax.step(th[t0:t1],rec[1,t0:t1],where='post',color='black',lw=0.9,label='λ$_{max}$(x)')
    v=rec[0]>rec[1]
    ax.fill_between(th[t0:t1],1500,2700,where=v[t0:t1],color='0.85',step='post',lw=0)
    ax.set_ylim(1500,2700); ax.set_ylabel('запит/с'); ax.set_title(tt,fontsize=9,loc='left')
axs[0].legend(loc='upper left',fontsize=8,frameon=False,ncol=2); axs[-1].set_xlabel('час доби, год')
fig.tight_layout(); fig.savefig('fig2.png')
