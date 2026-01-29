#!/usr/bin/env python3
"""
run_pipeline.py
===============
Master pipeline runner for TL1A–DcR3 cyclic peptide design.
Executes all four stages sequentially.
"""

import sys
import os

# Ensure we're in the project directory
os.chdir(os.path.dirname(os.path.abspath(__file__)))

print("╔" + "═" * 68 + "╗")
print("║  TL1A–DcR3 CYCLIC PEPTIDE DESIGN PIPELINE                        ║")
print("║  Based on PDB 3K51 crystal structure                              ║")
print("║  Target: Block TL1A–DcR3 PPI with NCAA-optimized cyclic peptide   ║")
print("╚" + "═" * 68 + "╝")

# Stage 1: Interface Analysis
print("\n\n" + "▶" * 3 + " STAGE 1: Interface Analysis")
from importlib import import_module
mod1 = __import__("01_interface_analysis")
summary = mod1.run_analysis()

# Stage 2: Cyclic Peptide Design
print("\n\n" + "▶" * 3 + " STAGE 2: Cyclic Peptide Design")
mod2 = __import__("02_cyclic_peptide_design")
designs = mod2.run_design()

# Stage 3: NCAA Optimization
print("\n\n" + "▶" * 3 + " STAGE 3: Non-Canonical Amino Acid Optimization")
mod3 = __import__("03_ncaa_optimization")
variants = mod3.run_optimization()

# Stage 4: Final Scoring & Ranking
print("\n\n" + "▶" * 3 + " STAGE 4: Final Scoring & Ranking")
mod4 = __import__("04_scoring_and_ranking")
candidates = mod4.run_final_scoring()

print("\n\n" + "═" * 70)
print("  PIPELINE COMPLETE")
print("═" * 70)
print("\nOutput files:")
print("  output/interface_analysis.json  — TL1A–DcR3 interface data")
print("  output/peptide_designs.json     — Cyclic peptide scaffolds")
print("  output/ncaa_optimized.json      — NCAA-optimized variants")
print("  output/final_candidates.json    — Final ranked candidates")
print()
