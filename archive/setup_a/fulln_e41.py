import sys, json, time, os
import pandas as pd
from bbcn import harness as H
from bbcn import repaired_branch as RB
RB.apply()
from bbcn import run_cde as CDE
COH={'TCGA':'data/binarized/tcga_brca_1082x135.csv',
     'METABRIC':'data/binarized/metabric_1980x135.csv',
     'ISPY2':'data/binarized/ispy2_988x135.csv'}
OUT='/home/claude/fulln_e41.json'
cohort=sys.argv[1]; a=int(sys.argv[2]); b=int(sys.argv[3])
df=pd.read_csv(COH[cohort],index_col=0); b=min(b,len(df))
acc=json.load(open(OUT)) if os.path.exists(OUT) else {}
d=acc.get(cohort, dict(cnt=0,apop=0,prolif=0,joint=0,dcommit=0,done=0))
t0=time.time()
for i in range(a,b):
    row={nd:int(df.iloc[i].get(nd,0)) for nd in df.columns}
    U=RB.patient_clamp_off(row); init=dict(row); init['CLAMP_OFF']=U
    init['PHLPP']=int(row.get('FOXO3',0) and not U); init['COMMIT']=0
    r=CDE.run_patient_cde(init,kernel_method='stabilize')
    ss=r['stage_summary']; ap=ss['Apoptosis_ON']['passed']
    d['cnt']+=1; d['apop']+=ap; d['prolif']+=ss['Proliferation_OFF']['passed']; d['joint']+=r['terminal_pass']
    if ap:
        fb=r['final_bus']
        d['dcommit']+=int(CDE._is_free_fixed_point(fb) and fb.get('CASP3',0)==1)
d['done']=max(d['done'],b)
acc[cohort]=d; json.dump(acc,open(OUT,'w'))
dt=time.time()-t0
print(f"{cohort} [{a}:{b}] {b-a} pts in {dt:.0f}s ({dt/max(1,b-a)*1000:.0f} ms/pt) | cohort total cnt={d['cnt']} done_to={d['done']}")
