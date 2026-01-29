#!/usr/bin/env python3
"""
02_cyclic_peptide_design.py
===========================
De novo design of cyclic peptides mimicking the DcR3 epitope that binds TL1A.

Design strategy:
  1. Extract the key epitope segments from DcR3 CRD2 and CRD3 that contact TL1A
  2. Design several cyclic peptide scaffolds:
     a) Single-loop: mimics the primary CRD2 hotspot loop (res 56-78)
     b) Bi-loop: links CRD2 and CRD3 hotspots via a short linker
     c) Miniprotein-like: includes β-hairpin mimicry from CRD2
  3. Add cyclization chemistry (head-to-tail, disulfide, thioether, CLIPS)
  4. Score designs by epitope coverage, predicted structure, and drug-likeness

Based on PDB 3K51 interface analysis from 01_interface_analysis.py.
"""

import json
import os
import itertools
from dataclasses import dataclass, field, asdict
from typing import List, Dict, Optional


# =============================================================================
# DATA FROM INTERFACE ANALYSIS
# =============================================================================

# Key DcR3 epitope sequences from 3K51 interface
# Segment 1: CRD2 loop (residues 56-92) — primary binding hotspot
CRD2_EPITOPE = "HCSPGQTCFPNHHKCDCYPGSYYGSSCQPACNAGTH"  # res 56-92
# Segment 2: CRD3 loop (residues 105-113) — secondary hotspot
CRD3_EPITOPE = "THDNQTFPG"  # res 105-113

# Hotspot residues from energy analysis (DcR3 numbering)
HOTSPOT_RESIDUES = {
    62: ("THR", -1.5),   # H-bond to TL1A E163
    67: ("HIS", -1.0),   # Aromatic with TL1A F168
    68: ("HIS", -1.5),   # H-bond to TL1A D169
    75: ("TYR", -1.5),   # H-bond to TL1A D205
    76: ("TYR", -1.5),   # H-bond to TL1A Q206
    78: ("SER", -1.5),   # H-bond to TL1A H208
    88: ("ASN", -1.5),   # H-bond to TL1A R223
    92: ("HIS", -1.0),   # Aromatic with TL1A H227
    107: ("ASP", -1.5),  # H-bond to TL1A N134
    108: ("ASN", -1.5),  # H-bond to TL1A I135
    110: ("THR", -1.5),  # H-bond to TL1A E137
    111: ("PHE", -0.5),  # Hydrophobic with TL1A E138
}


@dataclass
class CyclicPeptideDesign:
    """Represents a single cyclic peptide design candidate."""
    name: str
    sequence: str                    # Linear sequence before cyclization
    parent_segments: List[str]       # Which DcR3 segments it derives from
    cyclization: str                 # head-to-tail / disulfide / thioether / CLIPS
    linker: str = ""                 # Linker sequence if bi-loop
    length: int = 0
    hotspot_coverage: float = 0.0    # Fraction of hotspots covered
    mw_daltons: float = 0.0
    clogp_estimate: float = 0.0
    num_hbond_donors: int = 0
    num_hbond_acceptors: int = 0
    num_rotatable_bonds: int = 0
    oral_bioavail_score: float = 0.0  # 0-1 heuristic
    notes: str = ""

    def __post_init__(self):
        self.length = len(self.sequence)
        self._compute_properties()

    def _compute_properties(self):
        """Estimate molecular properties from sequence."""
        # Average amino acid MW ~110 Da; subtract 18 per peptide bond formed
        # Cyclization removes one additional water
        self.mw_daltons = round(self.length * 110 - (self.length) * 18, 1)

        # Approximate H-bond donors/acceptors from sequence
        polar = set("STYNQDHKRECW")
        self.num_hbond_donors = sum(1 for aa in self.sequence if aa in "STYNQHKRW")
        self.num_hbond_donors += self.length  # backbone NH (reduced by N-Me later)
        self.num_hbond_acceptors = sum(1 for aa in self.sequence if aa in "STYNQDEC")
        self.num_hbond_acceptors += self.length  # backbone CO

        # cLogP rough: hydrophobic residues contribute +, polar contribute -
        hydrophobic = set("AILMFWVP")
        self.clogp_estimate = round(
            sum(0.5 if aa in hydrophobic else -0.3 for aa in self.sequence), 2
        )

        # Rotatable bonds ~ 2 per residue (phi/psi) minus constraints from cyclization
        self.num_rotatable_bonds = max(0, self.length * 2 - 2)

        # Oral bioavailability heuristic (Veber + cyclic peptide rules)
        # MW < 800, HBD < 12 (or N-Me reduces), polar SA limited
        score = 1.0
        if self.mw_daltons > 1200:
            score -= 0.3
        elif self.mw_daltons > 800:
            score -= 0.1
        if self.num_hbond_donors > 15:
            score -= 0.3
        elif self.num_hbond_donors > 10:
            score -= 0.15
        if self.num_rotatable_bonds > 20:
            score -= 0.2
        # Cyclization bonus
        score += 0.1
        if self.length <= 12:
            score += 0.1  # Smaller = generally better permeability
        self.oral_bioavail_score = round(max(0.0, min(1.0, score)), 2)


# =============================================================================
# DESIGN FUNCTIONS
# =============================================================================

def design_crd2_single_loop():
    """
    Design 1: Single-loop cyclic peptide from CRD2 primary hotspot.
    Targets DcR3 residues 60-78 (the core TL1A-binding loop).
    This captures hotspots T62, H67, H68, Y75, Y76, S78.
    """
    # Core epitope: residues 60-78 of DcR3
    # GQTCFPNHHKCDCYPGSYYGSS
    # Trim to essentials and add cyclization-compatible termini
    designs = []

    # Design 1a: Minimal loop (14-mer), head-to-tail cyclization
    seq_1a = "QTCFPNHHKCDCYP"  # res 61-74, captures H67,H68 hotspots
    designs.append(CyclicPeptideDesign(
        name="CP-CRD2-min14",
        sequence=seq_1a,
        parent_segments=["CRD2_61-74"],
        cyclization="head-to-tail",
        hotspot_coverage=2 / len(HOTSPOT_RESIDUES),
        notes="Minimal CRD2 loop, captures H67/H68 hotspots"
    ))

    # Design 1b: Extended loop (19-mer), head-to-tail
    seq_1b = "GQTCFPNHHKCDCYPGSYY"  # res 60-78
    designs.append(CyclicPeptideDesign(
        name="CP-CRD2-ext19",
        sequence=seq_1b,
        parent_segments=["CRD2_60-78"],
        cyclization="head-to-tail",
        hotspot_coverage=6 / len(HOTSPOT_RESIDUES),
        notes="Extended CRD2 loop, captures T62,H67,H68,Y75,Y76,S78"
    ))

    # Design 1c: Disulfide-cyclized (replace terminal residues with Cys)
    seq_1c = "CQTCFPNHHKCDCYPGSYYC"  # Cys-Cys staple
    designs.append(CyclicPeptideDesign(
        name="CP-CRD2-SS20",
        sequence=seq_1c,
        parent_segments=["CRD2_60-78"],
        cyclization="disulfide",
        hotspot_coverage=6 / len(HOTSPOT_RESIDUES),
        notes="Disulfide-cyclized CRD2 loop, extra rigidity from 2 disulfides"
    ))

    return designs


def design_crd3_single_loop():
    """
    Design 2: Single-loop cyclic peptide from CRD3 secondary hotspot.
    Targets DcR3 residues 105-113 (THDNQTFPG).
    Captures hotspots D107, N108, T110, F111.
    """
    designs = []

    # Design 2a: CRD3 loop as-is, head-to-tail (9-mer — small, good permeability)
    seq_2a = "THDNQTFPG"
    designs.append(CyclicPeptideDesign(
        name="CP-CRD3-9mer",
        sequence=seq_2a,
        parent_segments=["CRD3_105-113"],
        cyclization="head-to-tail",
        hotspot_coverage=4 / len(HOTSPOT_RESIDUES),
        notes="CRD3 hotspot loop, very druggable size"
    ))

    # Design 2b: Extended with flanking residues for better fold
    seq_2b = "QHETHDNQTFPGQL"  # res 101-114
    designs.append(CyclicPeptideDesign(
        name="CP-CRD3-ext14",
        sequence=seq_2b,
        parent_segments=["CRD3_101-114"],
        cyclization="head-to-tail",
        hotspot_coverage=4 / len(HOTSPOT_RESIDUES),
        notes="Extended CRD3 loop with flanking for structure"
    ))

    return designs


def design_biloop_hybrid():
    """
    Design 3: Bi-loop hybrid linking CRD2 and CRD3 hotspot residues.
    Uses a short Gly-Pro or β-Ala linker to bridge the two epitope segments.
    This maximizes hotspot coverage.
    """
    designs = []

    # CRD2 core: hotspot residues T62, H67, H68 -> minimal seq around them
    crd2_core = "TCFPNHH"  # res 62-68
    # CRD3 core: hotspot D107, N108, T110, F111
    crd3_core = "DNQTF"    # res 107-111

    # Design 3a: Gly-Pro linked, head-to-tail (16-mer)
    seq_3a = crd2_core + "GP" + crd3_core + "GP"  # linkers at both ends
    designs.append(CyclicPeptideDesign(
        name="CP-Biloop-GP16",
        sequence=seq_3a,
        parent_segments=["CRD2_62-68", "CRD3_107-111"],
        cyclization="head-to-tail",
        linker="GP",
        hotspot_coverage=6 / len(HOTSPOT_RESIDUES),
        notes="Biloop hybrid, GP linkers, maximizes coverage"
    ))

    # Design 3b: Thioether-cyclized (ClAc-Cys chemistry, Suga-style)
    seq_3b = "C" + crd2_core + "GG" + crd3_core + "G"  # N-terminal ClAc-D-Cys
    designs.append(CyclicPeptideDesign(
        name="CP-Biloop-thioether17",
        sequence=seq_3b,
        parent_segments=["CRD2_62-68", "CRD3_107-111"],
        cyclization="thioether",
        linker="GG",
        hotspot_coverage=6 / len(HOTSPOT_RESIDUES),
        notes="Thioether macrocycle (RaPID-style), biloop hybrid"
    ))

    return designs


def design_beta_hairpin_mimic():
    """
    Design 4: β-hairpin mimic from CRD2 β-strand region (res 71-78).
    The CRD2 YPGSYYGS segment forms a β-strand that packs against TL1A.
    A β-hairpin peptide can mimic this extended conformation.
    """
    designs = []

    # β-strand: YPGSYYGS (res 71-78) + DPro-Gly turn
    # Mirror image β-hairpin
    strand1 = "YPGSYY"
    turn = "pG"  # D-Pro-Gly type II' turn (lowercase = D-amino acid in convention)
    strand2 = "NHHKCF"  # Anti-parallel from CRD2 loop (reversed segment ~67-62)

    seq_4a = strand1 + "PG" + strand2  # Will use D-Pro in NCAA optimization
    designs.append(CyclicPeptideDesign(
        name="CP-Hairpin14",
        sequence=seq_4a,
        parent_segments=["CRD2_71-78", "CRD2_62-67_rev"],
        cyclization="head-to-tail",
        hotspot_coverage=5 / len(HOTSPOT_RESIDUES),
        notes="β-hairpin mimic, DPro-Gly turn in NCAA step, captures Y75,Y76,H67,H68"
    ))

    return designs


def design_all_candidates():
    """Generate all cyclic peptide design candidates."""
    all_designs = []
    all_designs.extend(design_crd2_single_loop())
    all_designs.extend(design_crd3_single_loop())
    all_designs.extend(design_biloop_hybrid())
    all_designs.extend(design_beta_hairpin_mimic())
    return all_designs


# =============================================================================
# ANALYSIS AND REPORTING
# =============================================================================

def rank_designs(designs: List[CyclicPeptideDesign]):
    """
    Rank designs by a composite score:
      score = 0.4 * hotspot_coverage + 0.3 * oral_bioavail + 0.3 * size_score
    where size_score favors 8-16 residue range.
    """
    for d in designs:
        size_score = 1.0 if 8 <= d.length <= 16 else (0.7 if d.length <= 20 else 0.4)
        d._rank_score = (
            0.4 * d.hotspot_coverage
            + 0.3 * d.oral_bioavail_score
            + 0.3 * size_score
        )
    return sorted(designs, key=lambda d: d._rank_score, reverse=True)


def run_design():
    print("=" * 70)
    print("  CYCLIC PEPTIDE DESIGN — DcR3 EPITOPE MIMICS")
    print("=" * 70)

    designs = design_all_candidates()
    ranked = rank_designs(designs)

    print(f"\nGenerated {len(designs)} candidate designs:\n")
    print(f"{'Rank':<5} {'Name':<25} {'Len':<5} {'Cycl':<12} "
          f"{'MW':<8} {'HotCov':<8} {'OralBA':<8} {'Score':<8}")
    print("-" * 80)

    for i, d in enumerate(ranked, 1):
        print(f"{i:<5} {d.name:<25} {d.length:<5} {d.cyclization:<12} "
              f"{d.mw_daltons:<8.0f} {d.hotspot_coverage:<8.2f} "
              f"{d.oral_bioavail_score:<8.2f} {d._rank_score:<8.3f}")

    print("\n--- Top 3 Candidates for NCAA Optimization ---")
    for i, d in enumerate(ranked[:3], 1):
        print(f"\n  #{i} {d.name}")
        print(f"     Sequence:     {d.sequence}")
        print(f"     Length:       {d.length} aa")
        print(f"     Cyclization:  {d.cyclization}")
        print(f"     Segments:     {', '.join(d.parent_segments)}")
        print(f"     MW:           {d.mw_daltons:.0f} Da")
        print(f"     Hotspot cov:  {d.hotspot_coverage:.1%}")
        print(f"     Oral BA:      {d.oral_bioavail_score:.2f}")
        print(f"     Notes:        {d.notes}")

    # Save results
    os.makedirs("output", exist_ok=True)
    results = [asdict(d) for d in ranked]
    # Add rank scores
    for i, d in enumerate(ranked):
        results[i]["rank"] = i + 1
        results[i]["composite_score"] = round(d._rank_score, 4)
    with open("output/peptide_designs.json", "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n[✓] Designs saved to output/peptide_designs.json")

    return ranked


if __name__ == "__main__":
    run_design()
