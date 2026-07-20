#!/usr/bin/env python3
r"""
reproduce_all.py -- ONE command to rebuild the BBCN paper from code.

    python reproduce_all.py                 # rebuild preprint from committed numbers + figures (fast)
    python reproduce_all.py --figures regen # also rebuild every figure + the rules appendix
    python reproduce_all.py --numbers regen # also regenerate every number first (slow, full N)
    python reproduce_all.py --all           # full self-generating build (numbers + figures + paper)

What it does, in order:
  1. (optional, --numbers regen) regenerate every number -> results/numbers/{numbers.tex,
     numbers.json,all_numbers.csv} via run_chunked.run_all(), from the shipped CSVs.
  2. (optional, --figures regen) regenerate every figure + the rules appendix into preprint/
     via make_figures.py (Appendix A from pathways.py, the flowcharts, and the kernel figures).
  3. copy results/numbers/numbers.tex into preprint/ (the paper \input's it; paper == code).
  4. rebuild preprint/BBCN_preprint.pdf  (pdflatex, 3 passes) if a LaTeX engine is present.
  5. rebuild preprint/BBCN_preprint.docx (pandoc) if pandoc is present.
  Steps 3-4 degrade gracefully: if the tool is missing it prints an install hint and keeps
  the committed PDF/DOCX, so the command never hard-fails.

Expected numbers: results/LOCKED_NUMBERS.md   Universal table: results/numbers/all_numbers.csv
"""
import os, sys, shutil, subprocess, argparse

HERE = os.path.dirname(os.path.abspath(__file__))
PRE = os.path.join(HERE, "preprint")
NUM = os.path.join(HERE, "results", "numbers")

def _ensure(packages):
    missing = [pip for mod, pip in packages if not _has(mod)]
    if not missing: return
    print(f"[setup] installing {', '.join(missing)} ...")
    for args in (
        [sys.executable, "-m", "pip", "install", "--quiet", *missing],
        [sys.executable, "-m", "pip", "install", "--quiet", "--user", *missing],
        [sys.executable, "-m", "pip", "install", "--quiet", "--break-system-packages", *missing],
    ):
        try: subprocess.check_call(args); return
        except Exception: continue
    print("[setup] Could not auto-install. Please run:  pip install numpy pandas"); sys.exit(1)

def _has(mod):
    try: __import__(mod); return True
    except ImportError: return False

def banner(t): print("\n" + "=" * 74); print(t); print("=" * 74)

def have(tool): return shutil.which(tool) is not None

def regenerate_numbers():
    banner("STEP 1  regenerate every number (full cohorts) -> results/numbers/")
    sys.path.insert(0, HERE)
    import run_chunked
    run_chunked.run_all()      # processes all 27 tasks, then emits the three files

def regenerate_figures():
    banner("STEP 2  regenerate figures + rules appendix -> preprint/")
    r = subprocess.run([sys.executable, os.path.join(HERE, "make_figures.py")])
    if r.returncode != 0:
        print("[figures] generation reported an error; keeping committed figures.")

def sync_numbers_into_preprint():
    src = os.path.join(NUM, "numbers.tex")
    dst = os.path.join(PRE, "numbers.tex")
    if os.path.exists(src):
        shutil.copyfile(src, dst); print(f"[sync] {src} -> {dst}")
    else:
        print(f"[sync] {src} missing; using the numbers.tex already in preprint/")

def build_pdf():
    banner("STEP 3  rebuild PDF (pdflatex)")
    engine = "pdflatex" if have("pdflatex") else ("xelatex" if have("xelatex") else None)
    if not engine:
        print("[pdf] no LaTeX engine found. To rebuild: install TeX Live "
              "(e.g. `apt-get install texlive-latex-recommended`). Keeping committed PDF.")
        return
    for i in range(3):
        subprocess.run([engine, "-interaction=nonstopmode", "BBCN_preprint.tex"],
                       cwd=PRE, capture_output=True, text=True)
    print(f"[pdf] rebuilt preprint/BBCN_preprint.pdf with {engine}")

def build_docx():
    banner("STEP 4  rebuild DOCX (pandoc)")
    if not have("pandoc"):
        print("[docx] pandoc not found. To rebuild: install pandoc (https://pandoc.org). "
              "Keeping committed DOCX.")
        return
    r = subprocess.run(["pandoc", "BBCN_preprint.tex", "-o", "BBCN_preprint.docx",
                        "--resource-path=.:figures"], cwd=PRE, capture_output=True, text=True)
    print("[docx] rebuilt preprint/BBCN_preprint.docx" if r.returncode == 0
          else "[docx] pandoc error:\n" + r.stderr[-600:])

def setup_b_checks():
    banner("Setup B structural checks (no data needed)")
    B = os.path.join(HERE, "setup_b", "code"); sys.path.insert(0, B)
    import bbcn_switch as bs
    resting = dict(SRC=0, RHEB=1, IGF1R=1, RTK_up=1, CDKN2A=1, E2F1=0, ATM=0, ATR=0)
    stress  = dict(SRC=1, RHEB=1, IGF1R=1, RTK_up=1, CDKN2A=1, E2F1=1, ATM=1, ATR=1)
    surv = dict(AKT1=1, PTEN=0, MTOR=1, PDPK1=1, PIK3CA=1, MDM2=1, TP53=0, PHLPP=0, FOXO3=0)
    print("  resting, start SURVIVAL ->", bs.label(bs.simulate(surv, resting, CLAMP_OFF=0)),
          "(expect survival; only genotoxic stress flips it)")
    print("  stress,  start SURVIVAL ->", bs.label(bs.simulate(surv, stress, CLAMP_OFF=0)),
          "(expect apoptotic under sustained ATM/ATR/E2F1)")

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--numbers", choices=["committed", "regen"], default="committed",
                    help="committed = use results/numbers as shipped; regen = recompute (slow)")
    ap.add_argument("--figures", choices=["committed", "regen"], default="committed",
                    help="committed = use shipped figures; regen = rebuild figures + rules appendix")
    ap.add_argument("--all", action="store_true",
                    help="shortcut for --numbers regen --figures regen (full self-generating build)")
    a = ap.parse_args()
    if a.all:
        a.numbers, a.figures = "regen", "regen"
    pkgs = [("numpy", "numpy"), ("pandas", "pandas")]
    if a.figures == "regen":
        pkgs.append(("matplotlib", "matplotlib"))
    _ensure(pkgs)
    if a.numbers == "regen":
        regenerate_numbers()
    if a.figures == "regen":
        regenerate_figures()
    sync_numbers_into_preprint()
    build_pdf()
    build_docx()
    try: setup_b_checks()
    except Exception as e: print("  (Setup B checks skipped:", e, ")")
    banner("DONE")
    print("Preprint : preprint/BBCN_preprint.{tex,pdf,docx}")
    print("Numbers  : results/numbers/{numbers.tex,numbers.json,all_numbers.csv}")
    print("Expected : results/LOCKED_NUMBERS.md")
