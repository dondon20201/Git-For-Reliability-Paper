"""Draw Figures 1-3 of the paper from the committed results files.

Usage:
    python make_figures.py RESULTS_DIR [OUT_DIR]

RESULTS_DIR is the `expected/` folder of this repository, or its `output/`
folder after a full run. OUT_DIR defaults to the current folder. The script
writes fig1_reliability.pdf, fig2_projection.pdf and fig3_prior_studies.pdf.

Files read:
    results/step2/reliability_by_cell_construct.csv   (step 05)  Figures 1-3
    results/step4/reliability_extended.csv            (step 11)  Figures 1-3
    elevation/results/e0/variance_ceiling_final.csv   (step 19)  Figure 2

Figure 2 draws R(k) = s2_true / (s2_true + s2_err * k_ref / k) from the
committed variance components; the open circle marks k_0.7 from the same file.
The prior-study positions in Figure 3 (Walasek & Stewart, 2021, 64-100
gambles; Mkrtchian et al., 2023, about 120 design points, ICC .68 and .78) are
taken from data/prior_studies.csv in this repository.

Requires Python 3 with pandas, numpy and matplotlib.
"""
import os, sys
import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, NullLocator, FixedFormatter
if len(sys.argv) < 2:
    sys.exit(__doc__)
E = os.path.join(sys.argv[1], '')
OUT = os.path.join(sys.argv[2] if len(sys.argv) > 2 else '.', '')
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({'font.family':'serif','font.serif':['DejaVu Serif'],'font.size':9,'axes.labelsize':9.5,
 'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':0.7,'xtick.major.width':0.7,'ytick.major.width':0.7,
 'legend.frameon':False,'legend.fontsize':8.5,'pdf.fonttype':42})
ORANGE='#D55E00'; BLUE='#0072B2'; GREY='#8a8a8a'
def logx(ax,ticks,lo,hi):
    ax.set_xscale('log'); ax.set_xlim(lo,hi)
    ax.xaxis.set_major_locator(FixedLocator(ticks)); ax.xaxis.set_major_formatter(FixedFormatter([str(t) for t in ticks]))
    ax.xaxis.set_minor_locator(NullLocator())
def refs(ax,xr,labels=True):
    for y in (.5,.7):
        ax.axhline(y,color='#555',lw=0.6,ls=(0,(2,2)),zorder=0)
        if labels: ax.text(xr,y+0.012,f'r = {y:.1f}'.replace('0.','.'),fontsize=7.5,color='#555',ha='right',va='bottom')
s=pd.read_csv(E+'results/step2/reliability_by_cell_construct.csv'); s=s[s.axis=='i_exp']
x=pd.read_csv(E+'results/step4/reliability_extended.csv'); x=x[x.axis=='i_exp_ext']
def series(d,c):
    g=d[d.construct==c].dropna(subset=['r_sb']).sort_values('k'); return 2*g.k.values,g.r_sb.values,g.ci_lo_sb.values,g.ci_hi_sb.values
# ---------- Figure 1
fig,ax=plt.subplots(figsize=(6.3,4.2))
for d,ls,mk in ((s,'-','o'),(x,'--','^')):
    for c in ('b_ev','b_risk','b_loss'):
        k,r,_,_=series(d,c); ax.plot(k,r,ls=ls,color=GREY,lw=0.9,marker=mk,ms=2.5,zorder=1)
    for c,col in (('ev_sens_agn',BLUE),('loss_avoid_adj',ORANGE)):
        k,r,lo,hi=series(d,c)
        ax.fill_between(k,lo,hi,color=col,alpha=0.13,lw=0,zorder=2)
        ax.plot(k,r,ls=ls,color=col,lw=1.7,marker=mk,ms=3.6,zorder=3)
refs(ax,3.85)
for t in ax.texts: t.set_ha('left')
logx(ax,[4,6,8,10,14,20,26],3.7,28); ax.set_ylim(-0.08,0.75); ax.set_yticks([0,.2,.4,.6])
ax.set_xlabel('Total distinct problems (log scale)'); ax.set_ylabel('Split-half reliability (SB-corrected)')
from matplotlib.lines import Line2D
h=[Line2D([],[],color=ORANGE,lw=1.7,label='Loss avoidance (EV-adjusted)'),Line2D([],[],color=BLUE,lw=1.7,label='EV sensitivity'),
   Line2D([],[],color=GREY,lw=0.9,label='Model coefficients'),
   Line2D([],[],color='k',lw=1,ls='-',marker='o',ms=3,label='Strict pool'),Line2D([],[],color='k',lw=1,ls='--',marker='^',ms=3,label='Expanded pool')]
fig.legend(handles=h,loc='lower center',ncol=3,fontsize=8,columnspacing=1.4,handlelength=2.4,bbox_to_anchor=(0.5,0.0))
fig.tight_layout(rect=(0,0.12,1,1)); fig.savefig(OUT+'fig1_reliability.pdf'); plt.close()
# ---------- Figure 2
v=pd.read_csv(E+'elevation/results/e0/variance_ceiling_final.csv').set_index('construct')
fig,ax=plt.subplots(figsize=(6.3,3.7))
for c,col in (('loss_avoid_adj',ORANGE),('ev_sens_agn',BLUE)):
    st,se,kr,k7=v.loc[c,'sigma2_true'],v.loc[c,'sigma2_err'],v.loc[c,'k_ref'],v.loc[c,'k_0.7']
    R=lambda k: st/(st+se*kr/k)
    k1=np.linspace(4,26,200); ax.plot(k1,R(k1),color=col,lw=1.8)
    k2=np.linspace(26,k7,200); ax.plot(k2,R(k2),color=col,lw=1.8,ls='--')
    ax.plot([k7],[0.7],'o',mfc='white',mec=col,mew=1.5,ms=6,zorder=5)
    ax.annotate(f'{k7:.0f}',(k7,0.7),xytext=(-4,6),textcoords='offset points',ha='right',fontsize=8,color=col)
    for d,mk in ((s,'o'),(x,'^')):
        kk,r,_,_=series(d,c); ax.plot(kk,r,mk,color=col,ms=3.4,alpha=0.75,zorder=4)
ax.axhline(.7,color='#555',lw=0.6,ls=(0,(2,2)),zorder=0); ax.text(4.1,.712,'r = .7',fontsize=7.5,color='#555')
ax.axvline(26,color='#999',lw=0.6,zorder=0)
ax.text(25,0.03,'observed',ha='right',fontsize=7.5,color='#666'); ax.text(27,0.03,'projected',ha='left',fontsize=7.5,color='#666')
logx(ax,[4,6,10,14,26,40,60,84],3.7,100); ax.set_ylim(0,0.82)
ax.set_xlabel('Total distinct problems (log scale)'); ax.set_ylabel('Implied split-half reliability R(k)')
h=[Line2D([],[],color=ORANGE,lw=1.8,label='Loss avoidance (EV-adjusted)'),Line2D([],[],color=BLUE,lw=1.8,label='EV sensitivity (pooled across forms)'),
   Line2D([],[],color='k',lw=0,marker='o',ms=3.4,label='Observed, strict pool'),Line2D([],[],color='k',lw=0,marker='^',ms=3.4,label='Observed, expanded pool')]
ax.legend(handles=h,loc='lower right',ncol=1,bbox_to_anchor=(1.0,0.08))
fig.tight_layout(); fig.savefig(OUT+'fig2_projection.pdf'); plt.close()
# ---------- Figure 3
fig,ax=plt.subplots(figsize=(6.3,3.9))
ax.axvspan(4,26,color='#cfe8e2',alpha=0.55,lw=0,zorder=0)
ax.axvspan(64,100,color='#f3dfb3',alpha=0.7,lw=0,zorder=0)
ax.axvline(120,color='#444',lw=1,ls='--',zorder=1)
for d,ls,mk in ((s,'-','o'),(x,'--','^')):
    for c,col in (('ev_sens_agn',BLUE),('loss_avoid_adj',ORANGE)):
        k,r,_,_=series(d,c); ax.plot(k,r,ls=ls,color=col,lw=1.5,marker=mk,ms=3,zorder=3)
refs(ax,55)
ax.text(10.2,0.93,'CPC18 observed grid\n(4–26 problems)',ha='center',va='top',fontsize=8,color='#2b6b5f')
ax.text(80,0.93,'Walasek &\nStewart (2021)\n64–100 gambles',ha='center',va='top',fontsize=7.8)
ax.text(128,0.93,'Mkrtchian\net al. (2023)\n~120 points\nICC .68 / .78',ha='left',va='top',fontsize=7.8)
logx(ax,[4,6,10,20,26,64,120],3.7,260); ax.set_ylim(0,1.0)
ax.set_xlabel('Total problems / unique design points (log scale)'); ax.set_ylabel('Split-half reliability (SB-corrected)')
h=[Line2D([],[],color=ORANGE,lw=1.5,label='Loss avoidance'),Line2D([],[],color=BLUE,lw=1.5,label='EV sensitivity'),
   Line2D([],[],color='k',lw=1,marker='o',ms=3,label='Strict pool'),Line2D([],[],color='k',lw=1,ls='--',marker='^',ms=3,label='Expanded pool')]
fig.legend(handles=h,loc='lower center',ncol=4,fontsize=8,bbox_to_anchor=(0.5,0.0))
fig.tight_layout(rect=(0,0.08,1,1)); fig.savefig(OUT+'fig3_prior_studies.pdf'); plt.close()
print('Wrote fig1_reliability.pdf, fig2_projection.pdf and fig3_prior_studies.pdf to', OUT)
