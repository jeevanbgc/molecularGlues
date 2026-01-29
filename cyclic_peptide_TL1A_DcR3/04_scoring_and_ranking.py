#!/usr/bin/env python3
"""
04_scoring_and_ranking.py
=========================
Multi-objective scoring, final ranking, and synthesis recommendation for
NCAA-optimized cyclic peptide candidates targeting TL1A–DcR3 PPI.

Scoring dimensions:
  1. Binding affinity estimate (hotspot coverage + contact energy retention)
  2. Cell permeability prediction (MW, NMe count, HBD, PSA, cLogP)
  3. Proteolytic stability prediction (NCAA composition, backbone modifications)
  4. Synthetic accessibility (length, NCAA complexity, cyclization chemistry)
  5. Drug-likeness (beyond Rule of 5 for macrocycles)

Final output: ranked candidates with synthesis-ready specifications.
"""

import json
import os
from dataclasses import dataclass
from typing import List, Dict


# =============================================================================
# SCORING MODELS
# =============================================================================

# Molecular weight of each standard amino acid
AA_MW = {
    "G": 57.02, "A": 71.04, "V": 99.07, "L": 113.08, "I": 113.08,
    "P": 97.05, "F": 147.07, "W": 186.08, "M": 131.04, "S": 87.03,
    "T": 101.05, "C": 103.01, "Y": 163.06, "H": 137.06, "D": 115.03,
    "E": 129.04, "N": 114.04, "Q": 128.06, "K": 128.09, "R": 156.10,
}

# NCAA MW deltas
NCAA_MW_DELTA = {
    "NMe-L": 14, "NMe-G": 14, "NMe-A": 14, "NMe-F": 14,
    "NMe-T": 14, "NMe-Y": 14,
    "D-P": 0, "D-F": 0, "D-A": 0, "D-H": 0,
    "Aib": 14, "Cle": 12,
    "β3-hF": 14, "β3-hH": 14,
    "Nal": 40, "4F-Phe": 18, "Cha": 26,
    "Orn": -14, "Oic": 52,
}


def compute_exact_mw(sequence: str, substitutions: list) -> float:
    """Compute more accurate MW from sequence and substitutions."""
    # Sum individual AA MWs, subtract water for each peptide bond
    mw = sum(AA_MW.get(aa, 110.0) for aa in sequence)
    mw -= (len(sequence) - 1) * 18.015  # peptide bonds
    mw -= 18.015  # cyclization bond
    # Add NCAA MW deltas
    for sub in substitutions:
        mw += NCAA_MW_DELTA.get(sub.get("ncaa_code", ""), 0)
    return round(mw, 2)


def score_permeability_detailed(mw: float, num_nme: int, length: int,
                                 sequence: str) -> Dict[str, float]:
    """
    Detailed permeability scoring based on:
    - Molecular weight (< 700 ideal for oral; < 1000 for macrocycles)
    - N-methylation fraction (target 30-50% for cyclosporine-like permeability)
    - Hydrogen bond donor count (reduced by NMe and cyclization)
    - Estimated PSA (polar surface area)
    """
    scores = {}

    # MW score
    if mw <= 700:
        scores["mw_score"] = 1.0
    elif mw <= 1000:
        scores["mw_score"] = 1.0 - (mw - 700) / 600
    elif mw <= 1500:
        scores["mw_score"] = 0.5 - (mw - 1000) / 2000
    else:
        scores["mw_score"] = 0.1

    # N-methylation fraction score
    nme_frac = num_nme / length if length > 0 else 0
    if 0.25 <= nme_frac <= 0.50:
        scores["nme_score"] = 1.0
    elif 0.15 <= nme_frac < 0.25 or 0.50 < nme_frac <= 0.60:
        scores["nme_score"] = 0.7
    elif nme_frac > 0:
        scores["nme_score"] = 0.4
    else:
        scores["nme_score"] = 0.2

    # HBD count (backbone NH minus NMe sites, plus side-chain donors)
    polar_donors = sum(1 for aa in sequence if aa in "STYNQHKRW")
    backbone_hbd = length - num_nme  # each NMe removes one backbone NH
    total_hbd = polar_donors + backbone_hbd
    if total_hbd <= 6:
        scores["hbd_score"] = 1.0
    elif total_hbd <= 10:
        scores["hbd_score"] = 0.7
    elif total_hbd <= 15:
        scores["hbd_score"] = 0.4
    else:
        scores["hbd_score"] = 0.2

    # Composite permeability
    scores["permeability_total"] = round(
        0.35 * scores["mw_score"]
        + 0.35 * scores["nme_score"]
        + 0.30 * scores["hbd_score"], 3
    )
    return scores


def score_stability_detailed(num_d_aa: int, num_nme: int, num_other: int,
                              length: int) -> Dict[str, float]:
    """
    Proteolytic stability scoring.
    D-amino acids, N-methylation, and α,α-disubstituted AAs all reduce
    protease recognition.
    """
    scores = {}

    total_modified = num_d_aa + num_nme + num_other
    mod_fraction = total_modified / length if length > 0 else 0

    # Backbone modification coverage
    if mod_fraction >= 0.5:
        scores["modification_coverage"] = 1.0
    elif mod_fraction >= 0.3:
        scores["modification_coverage"] = 0.8
    elif mod_fraction >= 0.15:
        scores["modification_coverage"] = 0.5
    else:
        scores["modification_coverage"] = 0.2

    # D-amino acid bonus (especially effective)
    d_bonus = min(1.0, num_d_aa * 0.25)
    scores["d_aa_bonus"] = d_bonus

    # Cyclization inherent stability
    scores["cyclization_bonus"] = 0.3

    scores["stability_total"] = round(
        0.50 * scores["modification_coverage"]
        + 0.30 * scores["d_aa_bonus"]
        + 0.20 * scores["cyclization_bonus"], 3
    )
    return scores


def score_synthetic_accessibility(length: int, num_ncaa: int,
                                   cyclization: str) -> Dict[str, float]:
    """
    Synthetic accessibility scoring.
    Shorter peptides with fewer NCAAs and simpler cyclization are easier.
    """
    scores = {}

    # Length score
    if length <= 10:
        scores["length_score"] = 1.0
    elif length <= 14:
        scores["length_score"] = 0.8
    elif length <= 18:
        scores["length_score"] = 0.6
    else:
        scores["length_score"] = 0.4

    # NCAA complexity
    ncaa_frac = num_ncaa / length if length > 0 else 0
    if ncaa_frac <= 0.3:
        scores["ncaa_complexity"] = 1.0
    elif ncaa_frac <= 0.5:
        scores["ncaa_complexity"] = 0.7
    else:
        scores["ncaa_complexity"] = 0.4

    # Cyclization chemistry
    cycl_scores = {
        "head-to-tail": 0.9,
        "disulfide": 0.85,
        "thioether": 0.7,
        "CLIPS": 0.6,
        "lactam": 0.75,
    }
    scores["cyclization_ease"] = cycl_scores.get(cyclization, 0.5)

    scores["synthesis_total"] = round(
        0.40 * scores["length_score"]
        + 0.35 * scores["ncaa_complexity"]
        + 0.25 * scores["cyclization_ease"], 3
    )
    return scores


def score_drug_likeness(mw: float, hbd: int, hba: int, clogp: float,
                         rotatable_bonds: int) -> Dict[str, float]:
    """
    Beyond Rule of 5 (bRo5) drug-likeness for macrocycles.
    Based on Doak et al. (2016) analysis of oral macrocyclic drugs.
    """
    scores = {}

    # bRo5 MW: oral macrocycles up to ~1000 Da
    scores["bro5_mw"] = 1.0 if mw <= 1000 else max(0.0, 1.0 - (mw - 1000) / 500)

    # bRo5 HBD: < 12 for macrocycles
    scores["bro5_hbd"] = 1.0 if hbd <= 8 else max(0.0, 1.0 - (hbd - 8) / 8)

    # bRo5 cLogP: -2 to 6 range
    if -2 <= clogp <= 6:
        scores["bro5_clogp"] = 1.0
    else:
        scores["bro5_clogp"] = 0.5

    # Rotatable bonds: < 20 for macrocycles
    scores["bro5_rotbonds"] = 1.0 if rotatable_bonds <= 15 else max(
        0.0, 1.0 - (rotatable_bonds - 15) / 10)

    scores["druglikeness_total"] = round(
        0.30 * scores["bro5_mw"]
        + 0.30 * scores["bro5_hbd"]
        + 0.20 * scores["bro5_clogp"]
        + 0.20 * scores["bro5_rotbonds"], 3
    )
    return scores


# =============================================================================
# FINAL COMPOSITE SCORING
# =============================================================================

@dataclass
class FinalCandidate:
    rank: int
    name: str
    parent: str
    original_seq: str
    optimized_display: str
    cyclization: str
    mw: float
    length: int

    # Sub-scores
    affinity_score: float
    permeability_score: float
    stability_score: float
    synthesis_score: float
    druglikeness_score: float

    # Final
    final_score: float
    recommendation: str
    synthesis_spec: str


def generate_synthesis_spec(variant: dict) -> str:
    """Generate synthesis-ready specification string."""
    parts = []
    parts.append(f"Cyclization: {variant['cyclization']}")

    seq = list(variant["original_sequence"])
    ncaa_positions = {}
    for sub in variant.get("substitutions", []):
        ncaa_positions[sub["position"]] = sub["ncaa_code"]

    seq_parts = []
    for i, aa in enumerate(seq):
        if i in ncaa_positions:
            seq_parts.append(f"[{ncaa_positions[i]}]")
        else:
            seq_parts.append(aa)

    parts.append(f"Sequence: cyclo-({'–'.join(seq_parts)})")

    if variant["cyclization"] == "head-to-tail":
        parts.append("Chemistry: Fmoc SPPS, on-resin cyclization (PyBOP/DIPEA)")
    elif variant["cyclization"] == "disulfide":
        parts.append("Chemistry: Fmoc SPPS, oxidative folding (DMSO/buffer)")
    elif variant["cyclization"] == "thioether":
        parts.append("Chemistry: ClAc-D-aa at N-term, Cys thioether cyclization")

    return " | ".join(parts)


def run_final_scoring():
    print("=" * 70)
    print("  FINAL SCORING & RANKING — TL1A–DcR3 CYCLIC PEPTIDE CANDIDATES")
    print("=" * 70)

    # Load NCAA-optimized variants
    ncaa_file = "output/ncaa_optimized.json"
    if os.path.exists(ncaa_file):
        with open(ncaa_file) as f:
            variants = json.load(f)
    else:
        print("[!] Run 03_ncaa_optimization.py first.")
        return []

    candidates = []

    for v in variants:
        seq = v["original_sequence"]
        length = len(seq)
        mw = compute_exact_mw(seq, v.get("substitutions", []))
        num_nme = v.get("num_nme", 0)
        num_d_aa = v.get("num_d_aa", 0)
        num_other = v.get("num_other_ncaa", 0)

        # 1. Affinity (from NCAA optimization output)
        aff = v.get("affinity_retention", 0.5)

        # 2. Permeability
        perm = score_permeability_detailed(mw, num_nme, length, seq)

        # 3. Stability
        stab = score_stability_detailed(num_d_aa, num_nme, num_other, length)

        # 4. Synthetic accessibility
        total_ncaa = num_nme + num_d_aa + num_other
        synth = score_synthetic_accessibility(length, total_ncaa, v["cyclization"])

        # 5. Drug-likeness
        polar_donors = sum(1 for aa in seq if aa in "STYNQHKRW")
        hbd = polar_donors + length - num_nme
        hba = sum(1 for aa in seq if aa in "STYNQDEC") + length
        hydro = sum(0.5 if aa in "AILMFWVP" else -0.3 for aa in seq)
        rot_bonds = max(0, length * 2 - 2)
        dl = score_drug_likeness(mw, hbd, hba, hydro, rot_bonds)

        # Final composite
        final = round(
            0.30 * aff
            + 0.25 * perm["permeability_total"]
            + 0.20 * stab["stability_total"]
            + 0.10 * synth["synthesis_total"]
            + 0.15 * dl["druglikeness_total"], 3
        )

        # Recommendation
        if final >= 0.55:
            rec = "STRONG — prioritize for synthesis"
        elif final >= 0.45:
            rec = "MODERATE — synthesize in second tier"
        elif final >= 0.35:
            rec = "WEAK — consider with modifications"
        else:
            rec = "LOW — deprioritize"

        spec = generate_synthesis_spec(v)

        candidates.append(FinalCandidate(
            rank=0,
            name=v["variant_name"],
            parent=v["parent_name"],
            original_seq=seq,
            optimized_display=v.get("optimized_display", seq),
            cyclization=v["cyclization"],
            mw=mw,
            length=length,
            affinity_score=aff,
            permeability_score=perm["permeability_total"],
            stability_score=stab["stability_total"],
            synthesis_score=synth["synthesis_total"],
            druglikeness_score=dl["druglikeness_total"],
            final_score=final,
            recommendation=rec,
            synthesis_spec=spec,
        ))

    # Sort by final score
    candidates.sort(key=lambda c: c.final_score, reverse=True)
    for i, c in enumerate(candidates, 1):
        c.rank = i

    # Print results
    print(f"\n{'Rank':<5} {'Name':<30} {'MW':<7} {'Aff':<6} {'Perm':<6} "
          f"{'Stab':<6} {'Syn':<6} {'DL':<6} {'FINAL':<7} {'Rec':<30}")
    print("-" * 115)
    for c in candidates:
        print(f"{c.rank:<5} {c.name:<30} {c.mw:<7.0f} {c.affinity_score:<6.2f} "
              f"{c.permeability_score:<6.2f} {c.stability_score:<6.2f} "
              f"{c.synthesis_score:<6.2f} {c.druglikeness_score:<6.2f} "
              f"{c.final_score:<7.3f} {c.recommendation:<30}")

    # Detailed top candidates
    print(f"\n{'=' * 70}")
    print(f"  SYNTHESIS-READY CANDIDATES")
    print(f"{'=' * 70}")

    for c in candidates[:5]:
        print(f"\n  === Rank #{c.rank}: {c.name} ===")
        print(f"  Parent design:     {c.parent}")
        print(f"  Original sequence: {c.original_seq}")
        print(f"  Optimized:         {c.optimized_display}")
        print(f"  Length:            {c.length} residues")
        print(f"  MW:                {c.mw:.1f} Da")
        print(f"  Cyclization:       {c.cyclization}")
        print(f"  ")
        print(f"  Scores:")
        print(f"    Affinity retention:  {c.affinity_score:.3f}")
        print(f"    Permeability:        {c.permeability_score:.3f}")
        print(f"    Stability:           {c.stability_score:.3f}")
        print(f"    Synth accessibility: {c.synthesis_score:.3f}")
        print(f"    Drug-likeness:       {c.druglikeness_score:.3f}")
        print(f"    FINAL SCORE:         {c.final_score:.3f}")
        print(f"  ")
        print(f"  Recommendation: {c.recommendation}")
        print(f"  Synthesis spec: {c.synthesis_spec}")

    # Save final results
    os.makedirs("output", exist_ok=True)
    results = []
    for c in candidates:
        results.append({
            "rank": c.rank,
            "name": c.name,
            "parent": c.parent,
            "original_sequence": c.original_seq,
            "optimized_display": c.optimized_display,
            "cyclization": c.cyclization,
            "mw_daltons": c.mw,
            "length": c.length,
            "scores": {
                "affinity": c.affinity_score,
                "permeability": c.permeability_score,
                "stability": c.stability_score,
                "synthesis": c.synthesis_score,
                "drug_likeness": c.druglikeness_score,
                "FINAL": c.final_score,
            },
            "recommendation": c.recommendation,
            "synthesis_spec": c.synthesis_spec,
        })

    with open("output/final_candidates.json", "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n[✓] Final ranked candidates saved to output/final_candidates.json")

    return candidates


if __name__ == "__main__":
    run_final_scoring()
