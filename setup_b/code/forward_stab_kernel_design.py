import os
"""
Step 9: forward-stabilization kernel design (AIB's SR-paper method) applied
per-box to the three-box multirate switch.

Faithful to bbcn/stp.py: per box, build the box's L from rules, pick the box's
target sub-state z (the apoptotic projection onto that box), enumerate clamp-sets
S (size 1..2), project T@L, simulate forward with cycle detection, rank by
(impact -> delta -> steps -> causal). Kernels held per box's period.
"""
import numpy as np, gzip
from itertools import combinations
from collections import Counter
import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bbcn_switch import FAST, MID, SLOW, NODES, rules as switch_rules, simulate, label
import forward_stab as FS

# ---- the multirate boxes and their nodes ----
BOXES = {'FAST':FAST, 'MID':MID, 'SLOW':SLOW}
# druggable actuators (where clamps may be applied): FAST/MID signaling nodes only
DRUGGABLE_NODES = {'AKT1','PTEN','PIK3CA','PDPK1','MTOR'}  # FAST druggable actuators only (no direct MDM2 clamp)
# apoptotic target (full-state) -> project onto each box
APOP_TARGET = dict(AKT1=0,PTEN=1,MTOR=0,PDPK1=0,PIK3CA=0,MDM2=0,TP53=1,PHLPP=1,FOXO3=1)

def box_rule_fn(box, full_state, I, U):
    """Return a rule_fn(x_box, _) that updates ONLY the box's nodes, reading the
    rest of the state as frozen externals (the multirate 'other boxes held' view)."""
    def fn(x_box_list, _ext):
        s = dict(full_state)
        for nd,v in zip(box, x_box_list): s[nd]=v
        nxt = switch_rules(s, I, U)
        return [nxt[nd] for nd in box]
    return fn

def forward_select_box(box, full_state, I, U, method='stabilize'):
    """Forward-stabilization kernel on one box. Returns chosen clamp-set (node names) or None.

    method='ranked'    : preserved heuristic (impact->delta->steps->causal).
    method='stabilize' : algebraic Theorem-1 global-stabilization kernel (forward_stab.py),
                         over the same druggable candidate nodes, minimal cardinality.
    """
    bn = box
    x0 = [full_state[nd] for nd in bn]
    z  = [APOP_TARGET[nd] for nd in bn]
    if x0 == z: return ()   # already at target
    rule_fn = box_rule_fn(bn, full_state, I, U)
    n = len(bn)
    # candidate clamp indices: only druggable nodes in this box (0-based here, 1-based for FS)
    cand = [i for i,nd in enumerate(bn) if nd in DRUGGABLE_NODES]
    if not cand:
        return None
    if method == 'stabilize':
        rule_closed = lambda x: tuple(int(v) for v in rule_fn(list(x), None))
        sel = FS.stabilize_select_kernel(rule_closed, n, z,
                                         candidates=[i+1 for i in cand], max_k=len(cand))
        return tuple(bn[s-1] for s in sel) if sel else None
    # ---- ranked (preserved heuristic) ----
    results=[]
    for ksize in range(1, min(5,len(cand))+1):
        for Scomb in combinations(cand, ksize):
            # simulate forward: rules then pin S to target
            x=list(x0); matched=False; steps=0; visited=set()
            r=min(2**(n-len(Scomb)),200)
            for t in range(1,r+1):
                key=tuple(x)
                if key in visited: break
                visited.add(key)
                nxt=[int(v) for v in rule_fn(x,None)]
                for i in Scomb: nxt[i]=z[i]
                if nxt==z: matched=True; steps=t; break
                x=nxt
            if matched:
                impact=sum(1 for i in Scomb if x0[i]!=z[i])
                x1=[int(v) for v in rule_fn(x0,None)]
                for i in Scomb: x1[i]=z[i]
                delta=sum(1 for i in range(n) if x0[i]!=z[i] and x1[i]==z[i])
                causal=sum(i+1 for i in Scomb)
                results.append({'S':Scomb,'impact':impact,'delta':delta,'steps':steps,'causal':causal})
    if not results: return None
    # rank: max impact -> max delta -> min steps -> min causal
    mi=max(r['impact'] for r in results); R=[r for r in results if r['impact']==mi]
    md=max(r['delta'] for r in R); R=[r for r in R if r['delta']==md]
    ms=min(r['steps'] for r in R); R=[r for r in R if r['steps']==ms]
    mc=min(r['causal'] for r in R); R=[r for r in R if r['causal']==mc]
    chosen=R[0]['S']
    return tuple(bn[i] for i in chosen)

def design_patient_kernel(x0, I, U, method='stabilize'):
    """Run the patient to attractor; if resistant, design per-box forward-stab kernels.
    Returns (attractor, kernel_dict_by_box)."""
    final = simulate(x0, I, U)
    if label(final)=='APOPTOTIC':
        return 'APOPTOTIC', {}
    # design kernels per box, FAST first (fastest actuator), then MID
    kernels={}
    state = dict(final)   # design from the settled resistant state
    for bx in ['FAST','MID']:   # SLOW is target, not actuator
        sel = forward_select_box(BOXES[bx], state, I, U, method=method)
        if sel: kernels[bx]=sel
    # verify: apply the kernels (clamp) and re-simulate -> does it reach apoptosis?
    if kernels:
        clamp={nd:APOP_TARGET[nd] for bx in kernels for nd in kernels[bx]}
        xK=dict(final)
        for nd,v in clamp.items(): xK[nd]=v
        # held-clamp simulation
        s=dict(xK)
        for t in range(1,200+1):
            snap=dict(s); nf=switch_rules(snap,I,U)
            for nd in FAST: s[nd]=nf[nd]
            if t%5==0: s['MDM2']=switch_rules(snap,I,U)['MDM2']
            if t%25==0:
                for nd in SLOW: s[nd]=switch_rules(snap,I,U)[nd]
            for nd,v in clamp.items(): s[nd]=v   # HOLD the clamp every tick
        result=label(s)
    else:
        result=label(final)
    return result, kernels

# ---- minimal cohort harness (majority-gate Boolean init, §74) ----
SEED_SIG={'AKT1':['AKT1','MTOR','RPS6KB1','RPS6','EIF4EBP1','PDPK1','GSK3B'],'TP53':['CDKN1A','GADD45A','BBC3','MDM2','SESN1','RRM2B'],
 'MTOR':['MTOR','RPS6KB1','EIF4EBP1'],'PDPK1':['PDPK1','AKT1'],'PIK3CA':['PIK3CA','AKT1','PDPK1'],'PTEN':['PTEN'],'MDM2':['MDM2'],
 'PHLPP':['PHLPP1','PHLPP2','FOXO3'],'FOXO3':['FOXO3','FOXO1']}
AKT_ACT=['AKT1','MTOR','RPS6KB1','RPS6','EIF4EBP1']; FOXO_NUC=['FOXO3','FOXO1']
DAMAGE_SIG=['CHEK2','CHEK1','H2AFX','MDC1','RAD51','BRCA1','FANCD2','RNF168','TP53BP1','GADD45A','CDKN1A','ATM','ATR','MRE11','NBN','RAD50']
EXTRA=['SRC','RHEB','IGF1R','GRB2','IRS1','EGFR','ERBB2','CDKN2A','E2F1']
NEED=set(EXTRA+DAMAGE_SIG+AKT_ACT+FOXO_NUC)
for v in SEED_SIG.values(): NEED|=set(v)
def load_matrix(path, idcols=2, gz=False, header_idcols=None):
    """Load gene x sample matrix. header_idcols handles files whose header row
    has a different number of leading non-sample cells than data rows do (e.g.
    the I-SPY2 GEO matrix: 988 sample IDs in the header, but data rows are
    gene + 988 values). Pass header_idcols=0, idcols=1 for that file."""
    hc = idcols if header_idcols is None else header_idcols
    op=(lambda p:gzip.open(p,'rt')) if gz else open; rows={}
    with op(path) as f:
        samples=f.readline().rstrip('\n').split('\t')[hc:]
        for line in f:
            p=line.rstrip('\n').split('\t'); g=p[0].strip()
            if g in NEED:
                vals=[float(x) if x not in ('','NA','NaN') else np.nan for x in p[idcols:]]
                if len(vals)!=len(samples):
                    vals=(vals+[np.nan]*len(samples))[:len(samples)]
                rows[g]=np.array(vals)
    return samples, rows
def load_muts(path):
    tp53=set(); pik=set()
    with open(path) as f:
        line=f.readline()
        while line.startswith('#'): line=f.readline()
        h=line.rstrip('\n').split('\t'); gi=h.index('Hugo_Symbol'); bi=h.index('Tumor_Sample_Barcode'); vi=h.index('Variant_Classification')
        for line in f:
            p=line.rstrip('\n').split('\t')
            if len(p)<=max(gi,bi,vi): continue
            if p[gi]=='TP53' and p[vi] not in ('Silent',"3'UTR","5'UTR",'Intron','RNA'): tp53.add(p[bi])
            if p[gi]=='PIK3CA' and p[vi]=='Missense_Mutation': pik.add(p[bi])
    return tp53,pik
def binMedian(rows,n):
    B={}
    for g in NEED:
        if g in rows: v=rows[g]; med=np.nanmedian(v); b=(v>med).astype(int); b[np.isnan(v)]=0; B[g]=b
        else: B[g]=np.zeros(n,int)
    return B
def majority(B,genes,n):
    pres=[g for g in genes if g in B and B[g].any() or g in B]
    pres=[g for g in genes if g in B]
    if not pres: return np.zeros(n,int)
    cnt=np.sum([B[g] for g in pres],axis=0); thr=int(np.ceil(0.5*len(pres)))
    return (cnt>=thr).astype(int)

def run(name, samples, rows, tp53m=None, pikm=None, idtrim=lambda x:x, N=None, method='stabilize'):
    n=len(samples); B=binMedian(rows,n); N = N or n
    SEED={nd:majority(B,SEED_SIG[nd],n) for nd in SEED_SIG}
    akt=majority(B,AKT_ACT,n); fox=majority(B,FOXO_NUC,n); U=(akt&fox)
    damage=majority(B,DAMAGE_SIG,n)
    have=tp53m is not None
    routing=[]; ker_designed=0; ker_success=0; box_use=Counter(); node_use=Counter()
    for i in range(min(N,n)):
        sid=idtrim(samples[i]); dmg=int(damage[i])
        rtk=int(B['GRB2'][i] or B['IRS1'][i] or B.get('EGFR',np.zeros(n))[i] or B.get('ERBB2',np.zeros(n))[i])
        if have and sid in pikm: rtk=1
        I=dict(SRC=int(B['SRC'][i]),RHEB=int(B['RHEB'][i]),IGF1R=int(B['IGF1R'][i]),RTK_up=rtk,
               CDKN2A=int(B['CDKN2A'][i]),E2F1=int(B['E2F1'][i]),ATM=dmg,ATR=dmg)
        x0={nd:int(SEED[nd][i]) for nd in NODES}
        if have and sid in tp53m: x0['TP53']=0
        attr,kernels = design_patient_kernel(x0,I,int(U[i]),method=method)
        routing.append(attr)
        if label(simulate(x0,I,int(U[i])))=='SURVIVAL' or label(simulate(x0,I,int(U[i])))=='mixed':
            if kernels:
                ker_designed+=1
                if attr=='APOPTOTIC': ker_success+=1
                for bx,nds in kernels.items():
                    box_use[bx]+=1
                    for nd in nds: node_use[nd]+=1
    c=Counter(routing); m=min(N,n)
    flip=round(100*ker_success/max(ker_designed,1))
    print(f"\n{name} [method={method}]: N={m}")
    print(f"  post-design routing: APOP {100*c['APOPTOTIC']/m:.0f}% | SURV {100*c['SURVIVAL']/m:.0f}% | mixed {100*c['mixed']/m:.0f}%")
    print(f"  kernels designed for {ker_designed} non-apoptotic patients; flipped to APOP: {ker_success} ({flip}%)")
    print(f"  box usage: {dict(box_use)} | node usage: {dict(node_use.most_common())}")
    return dict(name=name, method=method, N=m, counts=dict(c),
                ker_designed=ker_designed, ker_success=ker_success, flip_pct=flip,
                node_use=dict(node_use.most_common()))

if __name__ == "__main__":
    # run small first to validate, then full
    tp53,pik=load_muts('data_mutations.txt'); s,r=load_matrix('data_mrna_seq_v2_rsem.txt')
    run('TCGA-BRCA (n=120 sample)',s,r,tp53,pik,idtrim=lambda x:x[:15],N=120)

    print("\n"+"="*70+"\nFULL COHORTS — forward-stabilization kernel design\n"+"="*70)
    tp53,pik=load_muts('data_mutations.txt'); s,r=load_matrix('data_mrna_seq_v2_rsem.txt')
    run('TCGA-BRCA',s,r,tp53,pik,idtrim=lambda x:x[:15])
    tp53m,pikm=load_muts('/tmp/brca_metabric/data_mutations.txt'); s,r=load_matrix('/tmp/brca_metabric/data_mrna_illumina_microarray.txt')
    run('METABRIC',s,r,tp53m,pikm)
    s,r=load_matrix('GSE194040_ISPY2ResID_AgilentGeneExp_990_FrshFrzn_meanCol_geneLevel_n988_txt.gz',idcols=1,gz=True,header_idcols=0)
    run('I-SPY2',s,r,None,None)
