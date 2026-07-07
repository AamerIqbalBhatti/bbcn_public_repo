#!/usr/bin/env python3
"""
make_figures.py -- regenerate every generated figure and the rules appendix into
the repository, from the shipped data, in one call. Each generator runs in its own
process for isolation (Setup A and Setup B ship same-named helper modules).

Outputs (all under preprint/):
  preprint/appendix_rules.tex                 (Appendix A, from setup_a/bbcn/pathways.py)
  preprint/figures/figC_flow_ranked.png       (Appendix B flowchart)
  preprint/figures/figC_flow_stabilize.png    (Appendix B flowchart)
  preprint/figures/figB_kernel_actuator.png   (Figure 4, Setup B kernel usage)
  preprint/figures/figA_kernel_usage.png      (Figure 5, Setup A kernel usage)
  preprint/figures/figA_pathway_kernel.png    (Figure 6, Setup A pathway x kernel)
"""
import os, sys, subprocess

HERE = os.path.dirname(os.path.abspath(__file__))


def ensure(mods):
    miss = []
    for mod, pip in mods:
        try:
            __import__(mod)
        except ImportError:
            miss.append(pip)
    if miss:
        print("[figures] installing", ", ".join(miss))
        for extra in ([], ["--user"], ["--break-system-packages"]):
            try:
                subprocess.check_call([sys.executable, "-m", "pip", "install", "--quiet", *extra, *miss])
                return
            except Exception:
                continue
        sys.exit("[figures] could not install: " + " ".join(miss))


def run(rel):
    print("\n[figures] >>>", rel)
    r = subprocess.run([sys.executable, os.path.join(HERE, rel)], cwd=HERE)
    if r.returncode != 0:
        sys.exit("[figures] FAILED: " + rel)


if __name__ == "__main__":
    ensure([("numpy", "numpy"), ("pandas", "pandas"), ("matplotlib", "matplotlib")])
    run("generate_rules_appendix.py")                       # Appendix A
    run("generate_flowcharts.py")                           # Appendix B
    run("make_setupB_kernel_heatmap.py")                    # Figure 4
    run(os.path.join("setup_a", "capture_setupA.py"))       # Setup A kernel capture (resumable)
    run(os.path.join("setup_a", "make_setupA_kernel_heatmaps.py"))  # Figures 5, 6
    run("make_tnbc_panel.py")                               # TNBC-like re-affirmation (Setup A)
    run("make_tnbc_setupB.py")                              # TNBC-like re-affirmation (Setup B)
    print("\n[figures] all figures and the rules appendix regenerated into preprint/")
