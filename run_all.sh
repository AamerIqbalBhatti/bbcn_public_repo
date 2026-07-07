#!/usr/bin/env bash
# Thin wrapper; the real entry point is run_all.py (portable, no bash needed).
cd "$(dirname "$0")" && exec python3 run_all.py "$@"
