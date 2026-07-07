"""
BBCN: Boolean Breast Cancer Network
============================================
Multi-Agent Boolean Output Regulation for SCC-Coupled Pathway Systems
BBCN Paper III — Aamer Iqbal Bhatti, KFUPM 2026

Package modules
---------------
pathways      : Boolean update rules for all 11 pathways (from getRules_v1.m)
kernels       : K_min (paper) and K_max (derived) kernel tables
derivation    : Exhaustive kernel derivation engine
simulation    : Patient-level Boolean dynamics simulator
phenotype     : Phenotype evaluation (5 stages + Terminal)
cohort        : TCGA cohort loader and preprocessor
survival      : Cox PH, Kaplan-Meier, log-rank analysis
plots         : All figure generators
"""

__version__ = "1.0.0"
__author__  = "Aamer Iqbal Bhatti"
__email__   = "aameriqbal.bhatti@kfupm.edu.sa"
