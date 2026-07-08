#!/usr/bin/env python3
"""Supplementary figure: the four computational pipelines of BBCN, what each one
takes in, the engine that runs it, the file it writes, and the paper element it
produces. Run: python3 fig_pipelines.py -> figures/fig_pipelines.png
"""
import os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

HERE=os.path.dirname(os.path.abspath(__file__)); os.makedirs(os.path.join(HERE,"figures"),exist_ok=True)
fig, ax = plt.subplots(figsize=(11.4, 7.2)); ax.set_xlim(0,100); ax.set_ylim(0,100); ax.axis("off")

def box(x,y,w,h,t,fc,ec="#555",fs=8.6,bold=False,tc="black",lw=1.1):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.3,rounding_size=1.2",fc=fc,ec=ec,lw=lw))
    ax.text(x+w/2,y+h/2,t,ha="center",va="center",fontsize=fs,fontweight=("bold" if bold else "normal"),color=tc,zorder=5)
def arr(x1,y1,x2,y2,c="#777",lw=1.5):
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle="-|>",mutation_scale=12,lw=lw,color=c,shrinkA=1,shrinkB=1))

# shared input
box(28,90,44,7,"Binarised cohorts\nTCGA 1082  \u00b7  METABRIC 1980  \u00b7  I-SPY2 988","#ededed",fs=9,bold=True)

# column headers
for x,t in [(6,"Pipeline (engine / key script)"),(52,"Writes"),(76,"Paper element")]:
    ax.text(x+ (19 if x>50 else 21),85.5,t,ha="center",fontsize=8.4,style="italic",color="#666")

rows=[
 ("#dbe7f3","#37628f","P1","Bare Setup A: three-tier whole-network\ncontroller on the 135-node network\n(generate_numbers.py)","results/numbers/\n(numbers.tex)","Table 1;\nSec 3.1\u20133.2;\nSupplement"),
 ("#dcefe0","#3c7a4e","P2","Repaired / delayed 136-node network:\nthree biological repairs + staged controller\n(generate_numbers2.py)","results/numbers2/\n(numbers2.tex)","Abstract;\nTable 3;\nSec 3.5"),
 ("#f7dede","#a23b3b","P3","Setup B: 9-node bistable switch on a\nmultirate schedule; per-patient routing\n(setup_b/code/cohort_pipeline.py)","routing +\nhysteresis","Table 2;\nFig 4\u20135;\nSec 3.3\u20133.4"),
 ("#faf0d7","#b8860b","P4","Switch vs whole-network durability:\nrelease-and-persist, minimal kernel\n(cde_vs_switch_cascade.py)","cde_vs_switch_\nsummary.csv\n(numbers_v2.tex)","Table 4;\nFig 7;\nSec 3.6"),
]
ys=[68,50,32,14]
for (fc,ec,tag,engine,out,paper),y in zip(rows,ys):
    ax.text(2.5,y+5.5,tag,ha="center",va="center",fontsize=10,fontweight="bold",color=ec)
    box(6,y,42,11,engine,fc,ec=ec,fs=8.2)
    box(52,y+1.5,20,8,out,"#f4f4f4",fs=8.0)
    box(76,y+1.5,20,8,paper,fc,ec=ec,fs=8.2,bold=True)
    arr(48,y+5.5,52,y+5.5); arr(72,y+5.5,76,y+5.5)
    arr(50,90,7,y+11,c="#bbb",lw=1.2)   # input bus to each engine

ax.text(50,4,"All four pipelines are required: each writes a different set of macros/figures that the manuscript reads.",
        ha="center",fontsize=8.2,style="italic",color="#555")
fig.tight_layout(pad=0.4)
out=os.path.join(HERE,"figures","fig_pipelines.png"); fig.savefig(out,dpi=200,bbox_inches="tight")
print("wrote",os.path.relpath(out,HERE))
