import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.colors as mc
from cips_ar import Params, AnalyticWindow
from retention import RetentionWindow
plt.rcParams.update({"font.family": "DejaVu Sans"})
fig = plt.figure(figsize=(8.85, 3.54), dpi=300)
fig.patch.set_facecolor("white")
fig.text(0.5, 0.955, "Floating-gate area ratio $A_R$ as a design knob for CuInP$_2$S$_6$ MFMIS transistors", ha="center", va="center", fontsize=10.5, fontweight="bold")
FACES=[]
def box(x0,x1,y0,y1,z0,z1,c):
    v=[(x0,y0,z0),(x1,y0,z0),(x1,y1,z0),(x0,y1,z0),(x0,y0,z1),(x1,y0,z1),(x1,y1,z1),(x0,y1,z1)]
    for q in [[0,1,2,3],[4,5,6,7],[0,1,5,4],[2,3,7,6],[1,2,6,5],[0,3,7,4]]: FACES.append(([v[k] for k in q], mc.to_rgba(c)))
def stack(ax, AR):
    FACES.clear(); s=np.sqrt(AR); h=0.5*s
    box(-1,-.5,-.4,1.4,0,.14,"#dddddd"); box(1.5,2,-.4,1.4,0,.14,"#dddddd"); box(-.5,1.5,-.4,1.4,0,.14,"#7fb2e5")
    box(0,1,0,1,.14,.34,"#a9d8b8"); box(0,1,0,1,.34,.42,"#e7c446")
    box(.5-h,.5+h,.5-h,.5+h,.42,.9,"#f4b4c4"); box(.5-h,.5+h,.5-h,.5+h,.9,.98,"#bbbbbb")
    ax.add_collection3d(Poly3DCollection([f for f,_ in FACES], facecolors=[c for _,c in FACES], edgecolors="k", linewidths=0.4))
    ax.set_xlim(-1,2); ax.set_ylim(-.6,1.6); ax.set_zlim(0,1.1); ax.set_box_aspect((3,2.2,1.1)); ax.view_init(32,-52); ax.set_axis_off()
for i,(a,lab) in enumerate([(0.25,"small $A_R$"),(2.0,"large $A_R$")]):
    ax=fig.add_axes([0.0, 0.50-0.46*i+0.02, 0.30, 0.42], projection="3d"); stack(ax,a)
    ax.set_title(lab+" = %g"%a, fontsize=8.5, y=0.88)
fig.text(0.15,0.03,"ferroelectric area = $A_R$ × channel area",ha="center",fontsize=7.5)
# middle: window vs AR
ax=fig.add_axes([0.385,0.20,0.28,0.62])
AR=np.logspace(-1.5,0.7,60)
for xp,c in [(6.0,"#0072B2")]:
    mw=[AnalyticWindow(Params(A_R=a,xp=xp)).MW() for a in AR]; ax.semilogx(AR,mw,color=c,lw=2,label="quasi-static, 6 V")
    ret=[RetentionWindow(Params(A_R=a,xp=xp),b=20).MW_t(3.15e8) for a in AR]; ax.semilogx(AR,ret,color="#D55E00",lw=2,label="10 yr (assumed $b$=20)")
ax.axvline(0.28,color="0.5",ls=":",lw=1); ax.text(0.31,0.85,"$A_R^{sat}$ ≈ 0.28",fontsize=7,color="0.3")
ax.text(0.035,1.42,"full window\n(plateau)",fontsize=7,va="top"); ax.text(2.7,1.25,"∝ 1/$A_R$",fontsize=7.5,color="#0072B2")
ax.set_xlabel("area ratio $A_R$",fontsize=8); ax.set_ylabel("memory window (V)",fontsize=8); ax.tick_params(labelsize=7)
ax.set_ylim(0,1.9); ax.legend(fontsize=6.5,frameon=False,loc="lower left")
for sp in ("top","right"): ax.spines[sp].set_visible(False)
# right: key results
t=["Key results (model study)","• Window saturates at its bound below\n   $A_R^{sat}$ ≈ 0.28 (6 V) and falls as $K/A_R$ above","• Drive needed to saturate scales ∝ $A_R$:\n   5.4 V at $A_R$=0.25 vs 18 V at $A_R$=1","• No interior optimum: traps and thermally\n   activated retention do not create one","• Published devices imply $A_R$ ≈ 0.2–0.5\n   (factor-of-two uncertainty); chip experiment\n   proposed to test"]
y=0.80
fig.text(0.70,0.84,t[0],fontsize=8.5,fontweight="bold",va="top")
for line in t[1:]:
    fig.text(0.70,y-0.05,line,fontsize=7.4,va="top"); y-=0.185 if "\n   " in line and line.count("\n")==2 else 0.155
fig.savefig("submission/graphical_abstract.png",dpi=300)
from PIL import Image
im=Image.open("submission/graphical_abstract.png").convert("RGB"); print(im.size); im.save("submission/graphical_abstract.tiff",dpi=(300,300),compression="tiff_lzw")
