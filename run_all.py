#!/usr/bin/env python3
"""
run_all.py -- one entry point to reproduce the BBCN results.

It auto-detects which cohorts are present in setup_a/data/binarized/ and runs the
four pipelines on whatever it finds. The public repository ships TCGA-BRCA only,
so a fresh clone reproduces the TCGA results; drop the METABRIC and I-SPY2
binarised inputs into the same folder (see DATA_ACCESS.md) and the identical code
reproduces all three. No code changes are needed either way.

Works on Windows, macOS and Linux with a plain Python install.

  python run_all.py                 # full run on every detected cohort
  python run_all.py --quick         # fast smoke test on a small sample (seconds)
  python run_all.py --pipelines 1,4 # run only selected pipelines
  python run_all.py --list          # show detected cohorts and exit

The four pipelines (each writes a different macro/figure the manuscript reads):
  1  Bare Setup A whole-network controller      -> results/numbers/  (numbers.tex)
  2  Repaired / delayed 136-node network        -> results/numbers2/ (numbers2.tex)   [SLOW]
  3  Setup B 9-node bistable switch             -> routing / hysteresis
  4  Switch vs whole-network durability         -> cde_vs_switch_summary.csv (numbers_v2.tex)
"""
import argparse
import glob
import os
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(ROOT, "setup_a", "data", "binarized")
COHORT_GLOBS = {
    "TCGA": "tcga_brca_*x135.csv",
    "METABRIC": "metabric_*x135.csv",
    "ISPY2": "ispy2_*x135.csv",
}

PIPELINES = {
    1: "Bare Setup A whole-network controller        (numbers.tex)",
    2: "Repaired / delayed 136-node network          (numbers2.tex)   [SLOW: full-cohort, minutes]",
    3: "Setup B 9-node bistable switch               (routing / hysteresis)",
    4: "Switch vs whole-network durability/minimality (numbers_v2.tex)",
}


def detect_cohorts():
    return [c for c, g in COHORT_GLOBS.items() if glob.glob(os.path.join(BIN, g))]


def switch_input_files():
    return sorted(glob.glob(os.path.join(ROOT, "setup_b", "data", "samples", "*_switch_inputs.tsv")))


def sh(args, cwd=ROOT):
    print("      $ python " + " ".join(args))
    sys.stdout.flush()
    return subprocess.call([sys.executable] + args, cwd=cwd)


def run_pipeline(n, quick):
    print("\n[{}/4] {}".format(n, PIPELINES[n]))
    if n == 1:
        return sh(["generate_numbers.py"])
    if n == 2:
        # reproduce_combined.py orchestrates the repaired runners; --quick samples N/cohort
        return sh(["reproduce_combined.py"] + (["--quick"] if quick else []))
    if n == 3:
        rc = 0
        files = switch_input_files()
        if not files:
            print("      (no *_switch_inputs.tsv present; skipping)")
            return 0
        for f in files:
            name = os.path.basename(f).split("_switch")[0].upper()
            rc |= sh(["setup_b/code/cohort_pipeline.py", "--expr", f, "--name", name])
        return rc
    if n == 4:
        return sh(["setup_a/cde_vs_switch_cascade.py"] + (["--sample", "40"] if quick else []))
    return 0


def main():
    ap = argparse.ArgumentParser(description="Reproduce BBCN results on whatever cohorts are present.")
    ap.add_argument("--quick", action="store_true", help="fast smoke test on a small sample")
    ap.add_argument("--pipelines", default="1,2,3,4", help="comma-separated subset, e.g. 1,4")
    ap.add_argument("--list", action="store_true", help="show detected cohorts and exit")
    args = ap.parse_args()

    cohorts = detect_cohorts()
    print("=" * 70)
    print("BBCN reproduction")
    print("Cohorts detected in setup_a/data/binarized/: " + (", ".join(cohorts) if cohorts else "NONE"))
    missing = [c for c in COHORT_GLOBS if c not in cohorts]
    if missing:
        print("Not present (controlled access; see DATA_ACCESS.md): " + ", ".join(missing))
    print("=" * 70)

    if args.list:
        return 0
    if not cohorts:
        print("No binarised cohort data found in setup_a/data/binarized/. Nothing to run.")
        return 1

    if args.quick:
        print("\nQuick smoke test (regression_check.py on a small sample)")
        print("Proves the pipelines are wired and reproduce the committed reference on the present cohort(s).")
        return sh(["regression_check.py", "--sample", "30"])

    try:
        selected = [int(x) for x in args.pipelines.split(",") if x.strip()]
    except ValueError:
        print("Bad --pipelines value; use e.g. --pipelines 1,4")
        return 2

    t0 = time.time()
    rc = 0
    for n in selected:
        if n not in PIPELINES:
            print("Unknown pipeline {}; valid are 1..4".format(n))
            continue
        rc |= run_pipeline(n, args.quick)

    print("\n" + "=" * 70)
    print("Done in {:.0f}s.  Numbers in results/ ; regenerate figures with make_figures.py.".format(time.time() - t0))
    print("Rebuild the manuscript:  cd paper3v2 && bash assemble.sh")
    print("=" * 70)
    return rc


if __name__ == "__main__":
    sys.exit(main())
