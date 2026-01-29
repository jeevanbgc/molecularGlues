# Cyclic Peptide Design: TL1A–DcR3 Interaction Blocker

## Overview

This project designs cyclic peptides that mimic the DcR3 (Decoy Receptor 3 / TNFRSF6B)
binding epitope to block the TL1A (TNFSF15)–DcR3 protein–protein interaction, based on the
crystal structure **PDB 3K51** (Liu et al., J. Biol. Chem. 2010).

## Biological Context

- **TL1A (TNFSF15)** is a TNF-superfamily cytokine that forms a homotrimer and signals
  through Death Receptor 3 (DR3/TNFRSF25), driving inflammation in IBD, rheumatoid
  arthritis, and other autoimmune diseases.
- **DcR3 (TNFRSF6B)** is a soluble decoy receptor that competitively blocks TL1A–DR3
  signalling. DcR3 is overexpressed in many cancers as an immune-evasion mechanism.
- A **cyclic peptide** mimicking DcR3's TL1A-binding epitope could serve as a defined,
  druggable inhibitor of TL1A signalling.

## Structure: PDB 3K51

| Property | Value |
|---|---|
| Title | Crystal structure of TL1A in complex with DcR3 |
| Resolution | 3.1 Å |
| Chains A,B,C | TL1A homotrimer (TNFSF15, residues 93-251) |
| Chains D,E,F | DcR3 (TNFRSF6B, CRD1-CRD4 domains) |
| Organism | Homo sapiens |
| Reference | Liu et al., J. Biol. Chem. 285, 15778-15786 (2010) |

## Pipeline

1. `01_interface_analysis.py` — Structural analysis of TL1A–DcR3 interface from 3K51
2. `02_cyclic_peptide_design.py` — De novo cyclic peptide design mimicking DcR3 epitope
3. `03_ncaa_optimization.py` — Non-canonical amino acid (NCAA) substitution for enhanced
   cell permeability and proteolytic stability
4. `04_scoring_and_ranking.py` — Multi-objective scoring and candidate ranking

## Non-Canonical Amino Acid Strategy

- **D-amino acids** at protease-susceptible sites
- **N-methylation** of backbone amides to improve membrane permeability
- **α,α-disubstituted amino acids (Aib)** for helix stabilization
- **β-amino acids** at select positions for protease resistance
- **Stapled residues** (olefin metathesis pairs) for conformational rigidity
