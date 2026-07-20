"""
repaired_branch.py — turn the REAL monolithic controller (run_cde.py) onto the
biologically-repaired network, by patching harness.PATHWAYS in place. Same three
repairs as honest_bbcn, now baked into the rules the controller reads:
  (1) p53-MDM2-ATM   (2) AKT1-FOXO3-PHLPP + CLAMP_OFF   (3) death-engaged commit latch + brake gating.
Adds 3 nodes: PHLPP, COMMIT, CLAMP_OFF. Original pathways.py untouched on disk.

EXECUTION PATTERN: callback installer. This module does NOT run the biology itself.
apply() builds rule functions (sg_rules, aux_rules, and the _override/_gate_brakes wrappers)
and STORES them into the pathway dicts (e.g. sg['rules'] = sg_rules). The harness/CDE spine
calls those stored rules later, once per node per step, when it runs a simulation. So this file
shows many 'def's and almost no direct calls by design: it INSTALLS the repair; the engine RUNS it.
The only thing run here is apply() itself, which other scripts call once to install the repair.
"""
from bbcn import harness as H

AKT_ACT=['AKT1','MTOR','RPS6KB1','EIF4EBP1']; FOXO_NUC=['FOXO3','FOXO1']
def patient_clamp_off(b):
    # CLAMP_OFF: failure of the AKT->14-3-3 nuclear-export clamp on FOXO3.
    # Normally AKT phosphorylates FOXO3 (T32/S253/S315), 14-3-3 binds and exports it from the
    # nucleus (Brunet 1999; Tzivion 2011, PMC3237389). The clamp fails under oxidative stress,
    # where stress kinases JNK and MST1 phosphorylate 14-3-3 and release FOXO3, so FOXO stays
    # nuclear despite active AKT (Sunayama/Gotoh 2005, PMC2171419; Sunters 2006 breast cancer).
    # The cohort cannot measure the cause (14-3-3 phospho-state), so CLAMP_OFF is defined by its
    # measurable EFFECT at baseline: strict majority of AKT nodes on AND strict majority of FOXO on.
    #   AKT  (4 nodes): at least 3 of 4 on.    FOXO (2 nodes): both on.  (Identical in honest_bbcn.)
    akt_on = sum(int(b.get(n,0)) for n in AKT_ACT)
    fox_on = sum(int(b.get(n,0)) for n in FOXO_NUC)
    return int(akt_on > len(AKT_ACT)/2 and fox_on > len(FOXO_NUC)/2)

def _override(pathway, overrides):
    orig=pathway['rules']; nodes=pathway['nodes']
    def wrapped(x,e):
        nxt=[int(v) for v in orig(x,e)]
        for nd,fn in overrides.items():
            if nd in nodes: nxt[nodes.index(nd)]=int(fn(x,e))
        return nxt
    pathway['rules']=wrapped

def _gate_brakes(pathway, brakes):
    orig=pathway['rules']; nodes=pathway['nodes']; idx=[nodes.index(b) for b in brakes if b in nodes]
    def wrapped(x,e):
        nxt=[int(v) for v in orig(x,e)]
        if int(e.get('COMMIT',0)):
            for i in idx: nxt[i]=0
        return nxt
    pathway['rules']=wrapped

def apply():
    P=H.PATHWAYS
    g=lambda e,k:int(e.get(k,0))
    # (1) p53 axis in AKT_Survival: MDM2, TP53 repaired; gate XIAP,CFLAR on COMMIT
    sv=P['AKT_Survival']  # nodes ['MDM2','TP53','JUN','XIAP','CFLAR']; x0=MDM2,x1=TP53
    _override(sv, {
        'MDM2': lambda x,e: (x[1] or g(e,'AKT1')) and not (g(e,'CDKN2A') and (g(e,'E2F1') or g(e,'ATM') or g(e,'ATR'))),
        'TP53': lambda x,e: (not x[0]) or g(e,'ATM') or g(e,'ATR'),
    })
    _gate_brakes(sv, ['XIAP','CFLAR'])
    # (2) AKT_Signaling: repair AKT1, FOXO3; add PHLPP as 7th node
    sg=P['AKT_Signaling']; nodes=sg['nodes']; orig=sg['rules']
    iA=nodes.index('AKT1'); iF=nodes.index('FOXO3'); sg['nodes']=nodes+['PHLPP']
    def sg_rules(x,e):
        nxt=[int(v) for v in orig(x[:len(nodes)],e)]
        phlpp_prev=x[len(nodes)] if len(x)>len(nodes) else 0
        foxo3_prev=x[iF]; akt_prev=x[iA]
        nxt[iA]=int((((g(e,'MTOR') and g(e,'PDPK1')) and (g(e,'PIK3CA') and not g(e,'PTEN'))) and not phlpp_prev) or (g(e,'CLAMP_OFF') and akt_prev))
        nxt[iF]=int((not akt_prev) or g(e,'CLAMP_OFF'))
        return nxt+[int(foxo3_prev and not g(e,'CLAMP_OFF'))]
    sg['rules']=sg_rules; sg['bus_exported']=sg.get('bus_exported',[])+['PHLPP']
    # (3) gate remaining brakes MCL1 (Apoptosis_Regulatory), BCL2/BCL2L1 (Apoptosis_Intrinsic)
    _gate_brakes(P['Apoptosis_Regulatory'], ['MCL1'])
    _gate_brakes(P['Apoptosis_Intrinsic'], ['BCL2','BCL2L1'])
    # new aux pathway: COMMIT latch (death-engaged) + U passthrough
    def aux_rules(x,e):
        return [int(x[0] or (g(e,'TP53') and not g(e,'AKT1'))), int(x[1])]
    P['RepairAux']={'nodes':['COMMIT','CLAMP_OFF'],'rules':aux_rules,'layer':'Downstream',
                    'role':'aux','sigma_order':99,'target':{},'bus_exported':['COMMIT','CLAMP_OFF']}
    # rebuild node tables
    H.ALL_NODES.clear() if hasattr(H.ALL_NODES,'clear') else None
    alln=[]; n2p={}
    for pn,pd in P.items():
        for nd in pd['nodes']:
            if nd not in n2p: alln.append(nd); n2p[nd]=pn
    H.ALL_NODES[:]=alln; H.NODE_TO_PATHWAY.clear(); H.NODE_TO_PATHWAY.update(n2p)
    # rules just changed to repaired -> invalidate externals + kernel caches so the
    # corrected, repaired-aware keys take effect (ledger E63).
    H._PATH_EXTERNALS.clear(); H._STAB_CACHE.clear()
    for _pd in P.values(): _pd.pop('externals', None)
    return len(alln)
