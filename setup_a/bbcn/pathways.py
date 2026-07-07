"""
bbcn.pathways  — ALL 24 pathways from getRules_v1.m (24 March 2026)
"""
from typing import List, Dict, Callable, Any
State = List[int]; Env = Dict[str, int]
def e_(e,k): return e.get(k,0)

def _rtk_egfr():
    def rules(x,e): return [
        int(x[0]),
        int((x[0] or (x[4] and x[5])) and not e_(e,'MAP2K1')),
        int(e_(e,'IRS1') or (e_(e,'LATS1') or (x[1] or e_(e,'SRC')))),
        int(x[3] and not e_(e,'PRKACA')),
        int(x[5] or x[0]),
        int(e_(e,'NRG1') or x[4]),
    ]
    return {'nodes':['EGF','EGFR','GRB2','GRB7','ERBB2','ERBB3'],'target':{},'role':'Upstream_RTK','bus_exported':['EGFR','GRB2','GRB7','ERBB2'],'rules':rules,'layer':'Upstream','sigma_order':1}

def _rtk_insulin():
    def rules(x,e): return [
        int(x[0]),
        int(x[0] or e_(e,'FOXO3')),
        int(not e_(e,'RPS6KB1') and (x[1] or e_(e,'IGF1R'))),
        int(e_(e,'GRB2') or x[2]),
        int(e_(e,'GRB2') or e_(e,'EGFR') or e_(e,'ERBB2') or x[2]),
    ]
    return {'nodes':['INS','INSR','IRS1','PIK3CA','ABL1'],'target':{},'role':'Upstream_RTK','bus_exported':['IRS1','PIK3CA','ABL1'],'rules':rules,'layer':'Upstream','sigma_order':2}

def _hormone():
    def rules(x,e): return [
        int(e_(e,'FOXO3') or not e_(e,'AKT1')),
        int(x[2]),
        int((e_(e,'FOXO3') or x[3]) and not e_(e,'AKT1')),
        int(e_(e,'PAX7')),
    ]
    return {'nodes':['AR','PGR','ESR1','KMT2D'],'target':{},'role':'Upstream_Hormone','bus_exported':['ESR1'],'rules':rules,'layer':'Upstream','sigma_order':3}

def _pi3k():
    def rules(x,e): return [
        int(not e_(e,'SRC')),
        int((e_(e,'PIK3CA') and not x[0]) or e_(e,'IGF1R')),
        int((e_(e,'PIK3CA') or e_(e,'DVL1')) and not e_(e,'AKT1')),
        int(e_(e,'HSP90AA1') or e_(e,'AKT1')),
    ]
    return {'nodes':['PTEN','PDPK1','RAC1','NOS3'],'target':{'PTEN':1},'role':'Terminal_support','bus_exported':['PTEN','PDPK1'],'rules':rules,'layer':'Intermediate','sigma_order':4}

def _mapk():
    def rules(x,e): return [
        int(x[4] and not x[5]),
        int(x[2]),
        int(x[0] or x[3]),
        int(x[4]),
        int(e_(e,'GRB2')),
        int((e_(e,'TP53') or e_(e,'FOXO3')) and not (x[7] or e_(e,'AKT1'))),
        int((e_(e,'FAS') or e_(e,'TRADD')) and not (e_(e,'AKT1') or x[7])),
        int(x[1] and not x[8]),
        int(x[7] or x[6] or e_(e,'JUN')),
    ]
    return {'nodes':['HRAS','MAP2K1','RAF1','KRAS','SOS1','NF1','MAPK8','MAPK1','DUSP1'],'target':{'MAPK1':0},'role':'Resistance','bus_exported':['MAP2K1','NF1','MAPK8','MAPK1','DUSP1'],'rules':rules,'layer':'Intermediate','sigma_order':5}

def _wnt():
    def rules(x,e): return [
        int(x[1]),
        int(not x[3] and e_(e,'GSK3B')),
        int(not (x[0] and x[1] and e_(e,'GSK3B')) and (x[7] or e_(e,'SRC') or e_(e,'YAP1'))),
        int(x[4] or e_(e,'NF1')),
        int(x[7] or e_(e,'CREB1')),
        int(e_(e,'GRB7') or x[3]),
        int(x[2] or e_(e,'JUN')),
        int(x[7]),
    ]
    return {'nodes':['APC','AXIN1','CTNNB1','DVL1','FZD3','RND1','TCF7L2','WNT1'],'target':{'CTNNB1':0},'role':'Invasion_OFF','bus_exported':['CTNNB1','DVL1','TCF7L2','WNT1'],'rules':rules,'layer':'Intermediate','sigma_order':6}

def _notch():
    def rules(x,e): return [
        int(e_(e,'ABL1')),
        int(x[1]),
        int(x[7]),
        int(x[7]),
        int(not e_(e,'MTOR')),
        int(not (x[2] or e_(e,'ABL1'))),
        int(not x[3] and x[5]),
        int(e_(e,'DLL1') or e_(e,'CTNNB1') or e_(e,'LCK') or e_(e,'MDM2')),
    ]
    return {'nodes':['ABL2','DLL1','HES1','HEY1','EIF4EBP1','MYOD1','MYOG','NOTCH1'],'target':{'NOTCH1':0},'role':'Invasion_OFF','bus_exported':['ABL2','EIF4EBP1'],'rules':rules,'layer':'Intermediate','sigma_order':7}

def _jak_stat():
    def rules(x,e): return [
        int(e_(e,'IL6') and e_(e,'IL6R')),
        int(x[0]),
        int(x[0] and not e_(e,'AKT1')),
        int(x[1] or x[2]),
    ]
    return {'nodes':['JAK2','STAT3','STAT1','SOCS3'],'target':{'STAT3':0},'role':'Xglobal_Tier1','bus_exported':['STAT3','STAT1'],'rules':rules,'layer':'Intermediate','sigma_order':8}

def _hippo():
    def rules(x,e): return [
        int(not x[1]),
        int(x[2] or e_(e,'LCK')),
        int(e_(e,'NEDD4L')),
    ]
    return {'nodes':['YAP1','LATS1','NF2'],'target':{'YAP1':0},'role':'Invasion_OFF','bus_exported':['YAP1','LATS1'],'rules':rules,'layer':'Intermediate','sigma_order':9}

def _akt_signaling():
    def rules(x,e): return [
        int((e_(e,'MTOR') and e_(e,'PDPK1')) and (e_(e,'PIK3CA') and not e_(e,'PTEN'))),
        int(not (x[2] or x[0])),
        int(e_(e,'MTOR') or e_(e,'PDPK1')),
        int(not (x[2] or e_(e,'PRKACA'))),
        int(not x[0]),
        int(not (x[0] or e_(e,'MAPK1'))),
    ]
    return {'nodes':['AKT1','FOXO1','SGK1','NEDD4L','FOXO3','GSK3B'],'target':{'AKT1':0,'FOXO3':1},'role':'Xglobal_Tier1','bus_exported':['AKT1','FOXO1','NEDD4L','FOXO3','GSK3B'],'rules':rules,'layer':'Intermediate','sigma_order':10}

def _mtor():
    def rules(x,e): return [
        int(not x[2]),
        int(x[0] or (e_(e,'PIK3CA') and not e_(e,'PTEN'))),
        int(not (e_(e,'AKT1') or e_(e,'MAPK1')) and (e_(e,'PRKAA2') or e_(e,'GSK3B'))),
        int(e_(e,'MAPK1') or e_(e,'AKT1') or e_(e,'MAPK8')),
    ]
    return {'nodes':['RHEB','MTOR','TSC2','TWIST1'],'target':{'TWIST1':0},'role':'Invasion_OFF','bus_exported':['MTOR'],'rules':rules,'layer':'Intermediate','sigma_order':11}

def _tf():
    def rules(x,e): return [
        int(e_(e,'ABL2') or x[1]),
        int(e_(e,'AKT1') and not e_(e,'PTEN')),
        int(not e_(e,'EIF4EBP1')),
        int(not e_(e,'STAT1')),
        int(not x[3]),
        int((e_(e,'CTNNB1') or e_(e,'PIM1') or e_(e,'PIM2') or e_(e,'PIM3') or
             e_(e,'MAPK1') or e_(e,'STAT3') or e_(e,'ESR1')) and not e_(e,'GSK3B')),
    ]
    return {'nodes':['CEBPB','CREB1','EIF4E','GATA3','CXCL8','MYC'],'target':{'MYC':0},'role':'Xglobal_Tier1','bus_exported':['CREB1','MYC'],'rules':rules,'layer':'Intermediate','sigma_order':12}

def _nfkb():
    def rules(x,e): return [
        int(e_(e,'TNF') or e_(e,'IL1A') or e_(e,'IL1B') or
            e_(e,'TLR2') or e_(e,'TLR4') or e_(e,'AKT1')),
        int(not x[0]),
        int(x[0] and not x[1]),
    ]
    return {'nodes':['IKBKB','NFKBIA','RELA'],'target':{'RELA':0},'role':'Resistance_OFF','bus_exported':['RELA'],'rules':rules,'layer':'Intermediate','sigma_order':13}

def _innate_immune():
    def rules(x,e): return [int(x[0]),int(x[1]),int(x[2]),int(x[3])]
    return {'nodes':['IL1A','IL1B','TLR2','TLR4'],'target':{},'role':'Inflammatory','bus_exported':['IL1A','IL1B','TLR2','TLR4'],'rules':rules,'layer':'Intermediate','sigma_order':14}

def _akt_survival():
    def rules(x,e): return [
        int(x[1] or x[0]),
        int(not x[0]),
        int((e_(e,'MAPK8') or e_(e,'MAPK1') or e_(e,'RELA')) and not (e_(e,'PPARG') or e_(e,'DUSP1'))),
        int(e_(e,'AKT1') or x[2]),
        int((e_(e,'AKT1') or e_(e,'RELA') or e_(e,'CREB1')) and x[1]),
    ]
    return {'nodes':['MDM2','TP53','JUN','XIAP','CFLAR'],'target':{'TP53':1,'XIAP':0,'CFLAR':0},'role':'Apoptosis_ON','bus_exported':['MDM2','TP53','JUN','XIAP','CFLAR'],'rules':rules,'layer':'Downstream','sigma_order':15}

def _apoptosis_regulatory():
    def rules(x,e): return [
        int(not e_(e,'GSK3B') and e_(e,'MAPK1')),
        int(e_(e,'TCF7L2')),
        int(e_(e,'AKT1') or e_(e,'WNT1')),
        int(x[2] or e_(e,'ABL1')),
    ]
    return {'nodes':['MCL1','PAK1','PRKACA','SRC'],'target':{'MCL1':0},'role':'Apoptosis_ON','bus_exported':['MCL1','PAK1','PRKACA','SRC'],'rules':rules,'layer':'Downstream','sigma_order':16}

def _apoptosis_extrinsic():
    def rules(x,e): return [
        int(x[1]),
        int(e_(e,'JUN') or e_(e,'RELA') or e_(e,'TP53')),
        int(x[0] or x[3]),
        int(x[0]),
        int(x[2] and not e_(e,'CFLAR')),
        int(x[4]),
        int((x[4] or x[7]) and not e_(e,'XIAP')),
        int(e_(e,'CYCS') and e_(e,'APAF1')),
    ]
    return {'nodes':['FAS','FASLG','FADD','TRADD','CASP8','BID','CASP3','CASP9'],'target':{'CASP3':1,'CASP9':1},'role':'Xglobal_Tier1','bus_exported':['FAS','TRADD'],'rules':rules,'layer':'Downstream','sigma_order':17}

def _apoptosis_intrinsic():
    def rules(x,e): return [
        int(not (e_(e,'PAK1') or e_(e,'AKT1') or e_(e,'MAPK1') or
                 e_(e,'PIM1') or e_(e,'PIM2') or e_(e,'PIM3'))),
        int((x[5] or e_(e,'TP53') or x[3]) and not (e_(e,'MCL1') or x[3] or x[4])),
        int((x[5] or e_(e,'TP53') or x[3]) and not (e_(e,'MCL1') or x[3] or x[4])),
        int(not (x[0] or x[5]) and (e_(e,'CREB1') or e_(e,'MAPK1'))),
        int(not (x[0] or x[5]) and e_(e,'CREB1')),
        int(e_(e,'FOXO3') or e_(e,'FOXO1')),
        int(x[1] or x[2]),
        int(e_(e,'TP53') or e_(e,'FOXO1') or e_(e,'FOXO3')),
    ]
    return {'nodes':['BAD','BAK1','BAX','BCL2','BCL2L1','BCL2L11','CYCS','APAF1'],'target':{'BAX':1,'CYCS':1,'APAF1':1},'role':'Apoptosis_ON','bus_exported':['BCL2','BCL2L1','CYCS','APAF1'],'rules':rules,'layer':'Downstream','sigma_order':18}

def _cellcycle():
    def rules(x,e): return [
        int(not (e_(e,'MYC') or e_(e,'SRC') or e_(e,'AKT1') or e_(e,'MAPK1'))),
        int((e_(e,'TCF7L2') or x[5] or e_(e,'MYC')) and not (x[0] or x[3] or e_(e,'GSK3B'))),
        int(not x[4]),
        int(e_(e,'TP53') or e_(e,'FOXO3') or not (e_(e,'MAPK1') or e_(e,'AKT1') or e_(e,'MYC'))),
        int(not (e_(e,'CDK4') or e_(e,'CDK6') or e_(e,'CDK2') or x[1])),
        int(e_(e,'YAP1')),
    ]
    return {'nodes':['CDKN2A','CCND1','E2F1','CDKN1A','RB1','TEAD1'],'target':{'CCND1':0,'E2F1':0},'role':'Proliferation_OFF','bus_exported':['CCND1'],'rules':rules,'layer':'Downstream','sigma_order':19}

def _cdk_cellcycle():
    def rules(x,e): return [
        int(e_(e,'CCND1') and x[4]),
        int(e_(e,'CCND1') and x[4]),
        int(x[0] or x[1]),
        int(x[0] or x[1]),
        int(e_(e,'FOXO1') or e_(e,'FOXO3')),
    ]
    return {'nodes':['CDK4','CDK6','CDK2','CCNE1','CDKN1B'],'target':{'CDK4':0,'CDK6':0},'role':'Proliferation_OFF','bus_exported':['CDK4','CDK6','CDK2'],'rules':rules,'layer':'Downstream','sigma_order':20}

def _dna_repair():
    def rules(x,e): return [int(x[0]),int(x[1]),int(x[0] or x[1]),int(x[2]),int(x[4]),int(x[0]),int(x[1]),int(x[2])]
    return {'nodes':['ATM','ATR','BRCA1','BRCA2','PARP1','CHEK1','CHEK2','RAD51'],'target':{},'role':'DNA_Repair','bus_exported':[],'rules':rules,'layer':'Downstream','sigma_order':21}

def _immune_checkpoint():
    def rules(x,e): return [
        int(e_(e,'LCK')),
        int(e_(e,'IFNG') or e_(e,'STAT3') or e_(e,'AKT1')),
    ]
    return {'nodes':['PDCD1','CD274'],'target':{},'role':'Immune_Checkpoint','bus_exported':[],'rules':rules,'layer':'Downstream','sigma_order':22}

def _aux_inputs1():
    def rules(x,e): return [int(x[i]) for i in range(8)]
    return {'nodes':['HSP90AA1','IGF1R','IL6','IL6R','LCK','NRG1','PAX7','IFNG'],'target':{},'role':'Boundary','bus_exported':['HSP90AA1','IGF1R','IL6','IL6R','LCK','NRG1','PAX7','IFNG'],'rules':rules,'layer':'Boundary','sigma_order':23}

def _aux_inputs2():
    def rules(x,e): return [int(x[i]) for i in range(7)]
    return {'nodes':['PIM1','PIM2','PIM3','PPARG','PRKAA2','RPS6KB1','TNF'],'target':{},'role':'Boundary','bus_exported':['PIM1','PIM2','PIM3','PPARG','PRKAA2','RPS6KB1','TNF'],'rules':rules,'layer':'Boundary','sigma_order':24}

PATHWAYS = {
    'RTK_EGFR':_rtk_egfr(),'RTK_Insulin':_rtk_insulin(),'Hormone':_hormone(),
    'PI3K':_pi3k(),'MAPK':_mapk(),'Wnt':_wnt(),'Notch':_notch(),
    'JAK_STAT':_jak_stat(),'Hippo':_hippo(),'AKT_Signaling':_akt_signaling(),
    'mTOR':_mtor(),'TF':_tf(),'NFkB':_nfkb(),'InnateImmune':_innate_immune(),
    'AKT_Survival':_akt_survival(),'Apoptosis_Regulatory':_apoptosis_regulatory(),
    'Apoptosis_Extrinsic':_apoptosis_extrinsic(),'Apoptosis_Intrinsic':_apoptosis_intrinsic(),
    'CellCycle':_cellcycle(),'CDK_CellCycle':_cdk_cellcycle(),
    'DNA_Repair':_dna_repair(),'Immune_Checkpoint':_immune_checkpoint(),
    'AuxInputs1':_aux_inputs1(),'AuxInputs2':_aux_inputs2(),
}

ACTIVE_PATHWAYS = {k:v for k,v in PATHWAYS.items() if v['layer']!='Boundary'}

SIGMA_BIO = [
    'RTK_EGFR','RTK_Insulin','Hormone','PI3K','MAPK','Wnt','Notch',
    'JAK_STAT','Hippo','AKT_Signaling','mTOR','TF','NFkB','InnateImmune',
    'AKT_Survival','Apoptosis_Regulatory','Apoptosis_Extrinsic',
    'Apoptosis_Intrinsic','CellCycle','CDK_CellCycle',
    'DNA_Repair','Immune_Checkpoint',
]

ALL_NODES = []
NODE_TO_PATHWAY = {}
for pname, pdef in PATHWAYS.items():
    for nd in pdef['nodes']:
        if nd not in NODE_TO_PATHWAY:
            ALL_NODES.append(nd)
            NODE_TO_PATHWAY[nd] = pname

def get_cancer_env():
    return {
        'MTOR':1,'PDPK1':1,'PIK3CA':1,'PTEN':0,'MAPK1':1,'PRKACA':0,
        'AKT1':1,'IL6':1,'IL6R':1,'MYC':1,'SRC':1,'GRB2':1,
        'TCF7L2':0,'GSK3B':0,'YAP1':0,'CDK4':1,'CDK6':1,'CDK2':1,
        'TP53':0,'FOXO3':0,'FOXO1':0,'MAPK8':0,'RELA':1,'PPARG':0,
        'DUSP1':0,'CYCS':0,'APAF1':0,'CFLAR':1,'XIAP':1,'BCL2':1,
        'BCL2L1':1,'MCL1':1,'ABL1':1,'DVL1':0,'IGF1R':0,'HSP90AA1':0,
        'WNT1':0,'STAT1':0,'STAT3':1,'EIF4EBP1':1,'ESR1':0,'PIM1':0,
        'PIM2':0,'PIM3':0,'CTNNB1':1,'ABL2':0,'JUN':0,'IL1A':0,'IL1B':0,
        'TNF':0,'CREB1':0,'NEDD4L':0,'BCL2L11':0,'BAD':0,'BAK1':0,
        'SOCS3':0,'CCND1':1,'E2F1':1,'MAP2K1':1,'RAF1':1,'NF1':0,
        'FAS':0,'TRADD':0,'LCK':0,'MDM2':1,'DLL1':0,'NOTCH1':0,
        'NRG1':0,'EGFR':1,'ERBB2':1,'IRS1':1,'GRB7':0,'LATS1':0,
        'PRKAA2':0,'RPS6KB1':0,'IFNG':0,'PAX7':0,'TLR2':0,'TLR4':0,
        'PAK1':0,'CCNE1':1,'NF2':0,'RAC1':0,'NOS3':0,'RHEB':1,
        'TSC2':0,'TWIST1':1,'HRAS':1,'KRAS':1,'SOS1':1,'NEDD4L':0,
        'SGK1':1,'FOXO1':0,
    }

def get_post_xglobal_env():
    e = get_cancer_env(); e.update({'AKT1':0,'MYC':0,'FOXO3':1,'CASP3':1}); return e

def get_healthy_env():
    e = {k:0 for k in get_cancer_env()}
    e.update({'PTEN':1,'TP53':1,'FOXO3':1,'GSK3B':1,'APAF1':1,
              'FOXO1':1,'AKT1':0,'MYC':0,'LATS1':1,'NF1':1,'DUSP1':1})
    return e

def get_all_envs():
    return [get_cancer_env(), get_post_xglobal_env(), get_healthy_env()]
