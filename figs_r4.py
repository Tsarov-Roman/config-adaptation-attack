import json,numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'Times New Roman','font.size':9,'axes.linewidth':0.6,'lines.linewidth':0.9,
  'xtick.direction':'in','ytick.direction':'in','savefig.dpi':300,'mathtext.fontset':'stix'})
O=json.load(open('results.json')); R=json.load(open('results_r4.json'))
fig,axs=plt.subplots(1,2,figsize=(6.3,2.6),sharey=True)
Ps=np.array([60,120,240,360,600,1200,1800,3600,7200])
for ax,T in zip(axs,[10,60]):
    for d,mk in [(0,'o'),(60,'s')]:
        ok=Ps/2>=2*T+2*d; pp=Ps[ok]
        ax.fill_between(pp,100*(T/2+d)/pp,100*(1.5*T+d)/pp,color='0.88' if d==0 else '0.78',alpha=0.7,lw=0)
        v1=[O[f'E3|T{T}|d{d}|P{P}|R1']['viol'][0] for P in Ps]
        ax.loglog(Ps,v1,mk+':',color='0.5',ms=3,mfc='white',label=f'R1, повна частка, d = {d} с')
        vi=[R['E7'][f'T{T}|d{d}|P{P}']['tinit'][0] for P in pp]
        ax.loglog(pp,vi,mk+'-',color='black',ms=3.5,label=f'R3, початкові порушення, d = {d} с')
    ax.set_title(f'T = {T} с',fontsize=9); ax.set_xlabel('період P, с'); ax.set_ylim(0.05,60)
axs[0].set_ylabel('частка часу, %'); axs[1].legend(fontsize=6.5,frameon=False,loc='lower left')
fig.tight_layout(); fig.savefig('fig3.png'); plt.close(fig)
# fig5
from trace import days
rs=np.load('rho_syn.npy')*100; rt=np.load('rho_tr.npy')*100
fig,axs=plt.subplots(1,2,figsize=(6.3,2.6))
for r,ls,lb in [(rs,'-','синтетичний профіль'),(rt,'--','траса WC98, 43 доби')]:
    xs=np.sort(r); axs[0].plot(xs,np.arange(len(xs))/len(xs),ls,color='black',label=lb)
axs[0].set_xlim(0,16); axs[0].set_xlabel('запас ρ$^{-}$, %'); axs[0].set_ylabel('частка доби'); axs[0].legend(fontsize=7.5,frameon=False,loc='lower right')
axs[0].set_title('а)',fontsize=9,loc='left')
bs=[-0.10,-0.05,-0.03,-0.02,-0.01,0.0]
for kind,ls in [('syn','-'),('tr','--')]:
    ex=[O[f'E6|{kind}|b{b}']['viol_exact'][0] for b in bs]; r3=[O[f'E6|{kind}|b{b}']['R3']['viol'][0] for b in bs]
    axs[1].plot([-100*b for b in bs],ex,'o'+ls,color='black',ms=3,label=('точна оцінка' if kind=='syn' else None))
    axs[1].plot([-100*b for b in bs],r3,'s'+ls,color='0.5',ms=3,mfc='white',label=('R3' if kind=='syn' else None))
axs[1].plot([],[],'-',color='black',label='синтетичний профіль'); axs[1].plot([],[],'--',color='black',label='траса WC98')
axs[1].set_xlabel('заниження оцінки |b|, %'); axs[1].set_ylabel('час порушення, %'); axs[1].legend(fontsize=7.5,frameon=False)
axs[1].set_title('б)',fontsize=9,loc='left')
fig.tight_layout(); fig.savefig('fig5.png'); plt.close(fig)
