import sys, time, json, os; sys.path.insert(0,'.')
import pandas as pd
from bbcn import harness as H
from bbcn import repaired_branch as RB
RB.apply()
from bbcn import run_cde as CDE
COH={'TCGA':'data/binarized/tcga_brca_1082x135.csv',
     'METABRIC':'data/binarized/metabric_1980x135.csv',
     'ISPY2':'data/binarized/ispy2_988x135.csv'}
name=sys.argv[1]; start=int(sys.argv[2]); end=int(sys.argv[3])
df=pd.read_csv(COH[name], index_col=0); n=len(df); end=n if end<0 else min(end,n)
t0=time.time(); apop=prolif=joint=dstrict=dcommit=0; cnt=0
for i in range(start,end):
    b={nd:int(df.iloc[i].get(nd,0)) for nd in df.columns}
    U=RB.patient_clamp_off(b)
    init=dict(b); init['CLAMP_OFF']=U; init['PHLPP']=int(b.get('FOXO3',0) and not U); init['COMMIT']=0
    r=CDE.run_patient_cde(init, kernel_method='stabilize'); cnt+=1
    ss=r['stage_summary']; a=ss['Apoptosis_ON']['passed']
    apop+=a; prolif+=ss['Proliferation_OFF']['passed']; joint+=r['terminal_pass']
    if a:
        fb=r['final_bus']; ffp=CDE._is_free_fixed_point(fb)
        dstrict+=int(ffp and fb.get('CASP3',0)==1 and fb.get('AKT1',1)==0)
        dcommit+=int(ffp and fb.get('CASP3',0)==1)
f=f'/home/claude/mono_{name}_{start}.json'
json.dump(dict(name=name,start=start,end=end,cnt=cnt,apop=apop,prolif=prolif,
               joint=joint,dstrict=dstrict,dcommit=dcommit,secs=round(time.time()-t0)),open(f,'w'))
print(json.load(open(f)))
