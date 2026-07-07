"""
bbcn.drugs
=============
Extended drug recommendation table with specific agents,
approval status, and mechanism of action.

Usage
-----
    from bbcn.drugs import DRUG_TABLE, get_drugs_for_pathway
    drugs = get_drugs_for_pathway('MAPK')
"""

DRUG_TABLE = [
    {"priority":1,"delta_H":3,"pathway":"Apoptosis_Intrinsic",
     "kernel":"BAX=1, CYCS=1, APAF1=1",
     "drug_class":"BCL-2 family inhibitor",
     "lead_drug":"Venetoclax (ABT-199)",
     "alternatives":["Navitoclax (ABT-263)","Obatoclax"],
     "approval":"FDA approved",
     "mechanism":"Inhibits BCL-2/BCL-XL, releasing BAX/BAK1 to trigger apoptosis"},

    {"priority":2,"delta_H":3,"pathway":"MAPK",
     "kernel":"MAP2K1=0, RAF1=0, MAPK1=0",
     "drug_class":"MEK + RAF inhibitor combination",
     "lead_drug":"Trametinib (MEK1/2)",
     "alternatives":["Dabrafenib (BRAF)","Cobimetinib (MEK)","Binimetinib"],
     "approval":"FDA approved (trametinib + dabrafenib)",
     "mechanism":"Blocks full MAPK cascade; prevents MYC re-activation"},

    {"priority":3,"delta_H":3,"pathway":"CellCycle",
     "kernel":"CCND1=0, E2F1=0, RB1=1",
     "drug_class":"CDK4/6 inhibitor",
     "lead_drug":"Palbociclib (Ibrance)",
     "alternatives":["Ribociclib (Kisqali)","Abemaciclib (Verzenio)"],
     "approval":"FDA approved (all three)",
     "mechanism":"Inhibits CDK4/6, stabilising RB1 to block E2F1 and CCND1"},

    {"priority":4,"delta_H":3,"pathway":"TF",
     "kernel":"MYC=0, CDKN2A=1, RB1=1",
     "drug_class":"BET bromodomain inhibitor",
     "lead_drug":"Molibresib (GSK525762)",
     "alternatives":["Birabresib (OTX015)","ZEN-3694","INCB054329"],
     "approval":"Phase I/II trials",
     "mechanism":"Suppresses MYC transcription by blocking BRD4 at MYC enhancers"},

    {"priority":5,"delta_H":3,"pathway":"AKT_Survival",
     "kernel":"TP53=1, XIAP=0, CFLAR=0",
     "drug_class":"MDM2 inhibitor + IAP antagonist",
     "lead_drug":"Idasanutlin (RG7388)",
     "alternatives":["Navtemadlin (AMG232)","Milademetan (DS-3032b)","Birinapant"],
     "approval":"Phase I/II trials",
     "mechanism":"MDM2i restores TP53=1; IAP antagonist removes XIAP/CFLAR block"},

    {"priority":6,"delta_H":2,"pathway":"AKT_Signaling",
     "kernel":"AKT1=0, FOXO3=1",
     "drug_class":"AKT inhibitor",
     "lead_drug":"Ipatasertib (GDC-0068)",
     "alternatives":["Capivasertib (AZD5363)","MK-2206","Uprosertib"],
     "approval":"Phase II/III",
     "mechanism":"Blocks AKT1 kinase; restores FOXO3 nuclear translocation"},

    {"priority":7,"delta_H":2,"pathway":"JAK_STAT",
     "kernel":"JAK2=0, STAT3=0",
     "drug_class":"JAK1/2 inhibitor",
     "lead_drug":"Ruxolitinib (Jakafi)",
     "alternatives":["Baricitinib","Fedratinib","Pacritinib"],
     "approval":"FDA approved (ruxolitinib for MF/PV)",
     "mechanism":"Inhibits JAK2 upstream, blocking STAT3 from IL6/IL6R"},

    {"priority":8,"delta_H":2,"pathway":"Apoptosis_Regulatory",
     "kernel":"MCL1=0, CASP9=1",
     "drug_class":"MCL1 inhibitor",
     "lead_drug":"S63845",
     "alternatives":["AZD5991","AMG-176","MIK665"],
     "approval":"Phase I trials",
     "mechanism":"Directly inhibits MCL1 anti-apoptotic protein"},

    {"priority":9,"delta_H":2,"pathway":"CDK_CellCycle",
     "kernel":"CDK4=0, CDK6=0",
     "drug_class":"CDK4/6 inhibitor (upstream)",
     "lead_drug":"Abemaciclib (Verzenio)",
     "alternatives":["Palbociclib (Ibrance)","Ribociclib (Kisqali)"],
     "approval":"FDA approved",
     "mechanism":"Addresses CDK4/6 upstream of CCND1/E2F1"},

    {"priority":10,"delta_H":2,"pathway":"Apoptosis_Extrinsic",
     "kernel":"CASP3=1, CASP9=1",
     "drug_class":"TRAIL receptor agonist",
     "lead_drug":"Dulanermin (rhTRAIL)",
     "alternatives":["Tigatuzumab (DR5 Ab)","Birinapant+TRAIL"],
     "approval":"Phase I/II trials",
     "mechanism":"Activates extrinsic apoptosis via CASP8→CASP3"},

    {"priority":11,"delta_H":1,"pathway":"PI3K",
     "kernel":"PTEN=1",
     "drug_class":"PI3K inhibitor",
     "lead_drug":"Alpelisib (BYL719)",
     "alternatives":["Copanlisib","Idelalisib","Buparlisib"],
     "approval":"FDA approved (alpelisib + fulvestrant)",
     "mechanism":"Inhibits PI3K to reduce PTEN suppression"},

    {"priority":12,"delta_H":1,"pathway":"NFkB",
     "kernel":"RELA=0",
     "drug_class":"NF-kB / proteasome inhibitor",
     "lead_drug":"Bortezomib (Velcade)",
     "alternatives":["Ixazomib (Ninlaro)","Carfilzomib (Kyprolis)"],
     "approval":"FDA approved (bortezomib for myeloma)",
     "mechanism":"Prevents IkB degradation, blocking RELA activation"},
]


def get_drugs_for_pathway(pathway_name: str) -> dict:
    """Return drug recommendation for a specific pathway."""
    return next((d for d in DRUG_TABLE
                 if d['pathway'] == pathway_name), None)


def get_approved_drugs() -> list:
    """Return only FDA-approved drug recommendations."""
    return [d for d in DRUG_TABLE if 'FDA' in d['approval']]


def format_patient_report(drug_recs: list) -> str:
    """Format drug recommendations as a readable text report."""
    lines = [
        "PERSONALISED DRUG RECOMMENDATIONS",
        "="*50,
        "Ranked by Hamming descent (ΔH) — higher = higher priority",
        "",
    ]
    for rec in drug_recs:
        pway  = rec['pathway']
        dH    = rec['delta_H']
        entry = get_drugs_for_pathway(pway)
        if entry is None:
            continue
        lines += [
            f"Priority #{rec['priority']}  ΔH={dH}",
            f"  Pathway    : {pway}",
            f"  Target     : {rec['kernel']}",
            f"  Drug class : {entry['drug_class']}",
            f"  Lead agent : {entry['lead_drug']}",
            f"  Approval   : {entry['approval']}",
            f"  Mechanism  : {entry['mechanism']}",
            "",
        ]
    return "\n".join(lines)


# ── Drug rationale generator ──────────────────────────────────
def generate_drug_rationale(report: dict, patient_row) -> str:
    """
    Full clinical drug rationale from algorithm output.
    Integrates pathway state, phenotype failures, AJCC stage,
    subtype, and patient condition into human-readable text.

    Parameters
    ----------
    report      : dict from simulation.single_patient_report()
    patient_row : pd.Series from cohort DataFrame

    Returns
    -------
    str : formatted rationale report
    """
    from bbcn.phenotype import check_all, vscc, PASS, FAIL, xglobal_violations

    pid    = report['patient_id']
    vscc_b = report['vscc_before']
    vscc_a = report['vscc_after']
    pheno  = report['phenotype_before']
    state  = report['init_state']

    subtype = patient_row.get('SUBTYPE', 'Unknown')
    stage   = patient_row.get('AJCC_PATHOLOGIC_TUMOR_STAGE', 'Unknown')
    age     = patient_row.get('AGE', '?')
    abl1    = state.get('ABL1', 0)
    xg_viols = xglobal_violations(state)

    lines = [
        "="*65,
        f"PERSONALISED DRUG RATIONALE — {pid}",
        "="*65,
        f"\nPATIENT: Subtype={subtype}  Stage={stage}  Age={age}",
        f"V_SCC:   {vscc_b} → {vscc_a}  (ΔV = {vscc_b-vscc_a})",
        f"ABL1:    {'=1 (basal-like, γp=0 for Apop_Reg)' if abl1 else '=0 (γp=1, monitor MCL1)'}",
        "",
        "PHENOTYPE BEFORE TREATMENT:",
    ]
    for s, v in pheno.items():
        lines.append(f"  {'✓' if v==PASS else '✗'} {s}: {v}")

    if xg_viols:
        lines += ["", "XGLOBAL VIOLATIONS (addressed first):"]
        for nd, cv in xg_viols.items():
            tgt = 0 if nd in ('AKT1','MYC') else 1
            lines.append(f"  {nd}: {cv} → {tgt}")

    lines += ["", "DRUG RECOMMENDATIONS WITH RATIONALE:", "─"*65]

    for rec in report['drug_priority']:
        if rec['delta_H'] == 0:
            continue
        entry = get_drugs_for_pathway(rec['pathway'])
        if not entry:
            continue

        lines += [
            f"\n#{rec['priority']} {rec['pathway']}  ΔH={rec['delta_H']}",
            f"  Kernel  : {rec['kernel']}",
            f"  Drug    : {entry['lead_drug']}",
            f"  Also    : {', '.join(entry['alternatives'][:2])}",
            f"  Approval: {entry['approval']}",
            f"  Mechanism: {entry['mechanism']}",
        ]

        # Stage-specific note
        if 'IV' in str(stage):
            lines.append(f"  Stage IV note: prefer CNS-penetrant agents where available")
        if 'Basal' in str(subtype) and rec['pathway'] == 'AKT_Signaling':
            lines.append(f"  Subtype note: AKT hyperactivated in Basal — strong indication")

    lines += ["", "="*65,
              f"Administer in σ_bio order. Skip ΔH=0 steps (already at target).",
              "="*65]

    return "\n".join(lines)
