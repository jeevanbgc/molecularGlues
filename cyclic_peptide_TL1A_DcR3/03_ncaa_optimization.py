#!/usr/bin/env python3
"""
03_ncaa_optimization.py
=======================
Non-canonical amino acid (NCAA) optimization of cyclic peptide candidates
to improve:
  - Cell membrane permeability (N-methylation, lipophilic NCAAs)
  - Proteolytic stability (D-amino acids, β-amino acids, Aib)
  - Conformational rigidity (stapled residues, Aib, Pro surrogates)
  - Target affinity (halogenated aromatics, extended side chains)

NCAA Library used:
  D-aa:    D-isomers of natural amino acids (protease resistance)
  NMe-aa:  N-methylated backbone (permeability, reduce HBD count)
  Aib:     α-aminoisobutyric acid (helix/turn stabilization)
  β3-aa:   β3-homo amino acids (protease resistance, backbone extension)
  Nal:     β-(2-naphthyl)-L-alanine (enhanced aromatic contacts)
  Cle:     cycloleucine (α,α-disubstituted, proteolysis resistant)
  hArg:    homoarginine (extended cation for salt bridges)
  Orn:     ornithine (shorter Lys analog)
  Cha:     cyclohexylalanine (hydrophobic, metabolically stable)
  4-F-Phe: 4-fluorophenylalanine (tuned aromatic electronics)
  Oic:     octahydroindole-2-carboxylic acid (rigid Pro replacement)

Strategy follows principles from:
  - Cyclosporine A (N-methylation pattern for oral bioavailability)
  - Sanguinamide B (alternating NMe pattern)
  - RaPID system (thioether macrocycles with NCAAs)
"""

import json
import os
from dataclasses import dataclass, field, asdict
from typing import List, Dict, Tuple, Optional
from copy import deepcopy


# =============================================================================
# NCAA DEFINITIONS
# =============================================================================

@dataclass
class NonCanonicalAA:
    """Definition of a non-canonical amino acid."""
    code: str           # Short code (e.g., "NMe-L", "D-H", "Aib")
    full_name: str
    category: str       # D-amino, N-methyl, alpha-disubstituted, beta-amino, modified
    replaces: str       # Which natural AA(s) it can replace ("*" = any)
    mw_delta: float     # MW change vs natural AA (Da)
    permeability_bonus: float   # 0-1 improvement to permeability
    stability_bonus: float      # 0-1 improvement to proteolytic stability
    affinity_effect: float      # -1 to +1 (neg = reduces, pos = enhances)
    notes: str = ""

NCAA_LIBRARY = [
    # --- D-amino acids ---
    NonCanonicalAA("D-P", "D-Proline", "D-amino", "P", 0, 0.0, 0.6, -0.1,
                   "Induces type II' β-turn; essential for hairpin designs"),
    NonCanonicalAA("D-F", "D-Phenylalanine", "D-amino", "F", 0, 0.05, 0.6, -0.2,
                   "Protease-resistant aromatic"),
    NonCanonicalAA("D-A", "D-Alanine", "D-amino", "A", 0, 0.0, 0.5, -0.1,
                   "Small, flexible D-residue"),
    NonCanonicalAA("D-H", "D-Histidine", "D-amino", "H", 0, 0.0, 0.6, -0.15,
                   "Protease-resistant His replacement"),

    # --- N-methylated amino acids ---
    NonCanonicalAA("NMe-L", "N-methyl-Leucine", "N-methyl", "L", 14, 0.3, 0.4, 0.0,
                   "Key for cyclosporine-like permeability"),
    NonCanonicalAA("NMe-G", "N-methyl-Glycine (sarcosine)", "N-methyl", "G", 14, 0.3, 0.3, 0.0,
                   "Sarcosine — excellent permeability enhancer"),
    NonCanonicalAA("NMe-A", "N-methyl-Alanine", "N-methyl", "A", 14, 0.3, 0.4, 0.0,
                   "Small N-methyl residue"),
    NonCanonicalAA("NMe-F", "N-methyl-Phenylalanine", "N-methyl", "F", 14, 0.3, 0.4, 0.05,
                   "N-methyl aromatic, maintains contact"),
    NonCanonicalAA("NMe-T", "N-methyl-Threonine", "N-methyl", "T", 14, 0.25, 0.3, -0.05,
                   "N-methylated polar — removes one HBD"),
    NonCanonicalAA("NMe-Y", "N-methyl-Tyrosine", "N-methyl", "Y", 14, 0.25, 0.35, -0.05,
                   "N-methylated Tyr"),

    # --- α,α-disubstituted ---
    NonCanonicalAA("Aib", "α-aminoisobutyric acid", "alpha-disubstituted", "A", 14, 0.1, 0.7, 0.0,
                   "Strong 3₁₀/α-helix inducer, protease resistant"),
    NonCanonicalAA("Cle", "Cycloleucine", "alpha-disubstituted", "L", 12, 0.1, 0.7, 0.0,
                   "Cyclopentane ring on Cα, very protease resistant"),

    # --- β-amino acids ---
    NonCanonicalAA("β3-hF", "β3-homo-Phenylalanine", "beta-amino", "F", 14, 0.05, 0.8, -0.1,
                   "Extra backbone methylene, highly protease resistant"),
    NonCanonicalAA("β3-hH", "β3-homo-Histidine", "beta-amino", "H", 14, 0.0, 0.8, -0.15,
                   "β-His analog for hotspot positions"),

    # --- Modified aromatics ---
    NonCanonicalAA("Nal", "β-(2-naphthyl)-L-alanine", "modified", "F,Y,W", 40, 0.1, 0.3, 0.2,
                   "Larger aromatic surface, enhanced π-stacking with TL1A"),
    NonCanonicalAA("4F-Phe", "4-fluoro-L-phenylalanine", "modified", "F", 18, 0.15, 0.2, 0.15,
                   "Fluorine enhances metabolic stability and aromatic interactions"),
    NonCanonicalAA("Cha", "cyclohexylalanine", "modified", "L,I,V", 26, 0.15, 0.4, 0.1,
                   "Metabolically stable hydrophobic replacement"),

    # --- Modified basic/polar ---
    NonCanonicalAA("Orn", "Ornithine", "modified", "K", -14, 0.0, 0.2, -0.05,
                   "Shorter Lys — still cationic, less flexible"),
    NonCanonicalAA("Oic", "Octahydroindole-2-carboxylic acid", "modified", "P", 52, 0.15, 0.5, 0.1,
                   "Rigid bicyclic Pro replacement, constrains backbone"),
]


# =============================================================================
# SUBSTITUTION RULES
# =============================================================================

# Position classification for NCAA substitution strategy
POSITION_ROLES = {
    "hotspot": "Residue at binding interface hotspot — preserve side chain",
    "scaffold": "Structural residue for peptide fold — can modify backbone",
    "solvent": "Solvent-exposed, non-contact — best site for NMe/D-aa",
    "turn": "Turn/linker residue — ideal for D-Pro/Aib",
}


@dataclass
class NCASubstitution:
    """A single NCAA substitution at a specific position."""
    position: int           # 0-indexed in peptide sequence
    original_aa: str        # Single letter
    ncaa: NonCanonicalAA
    role: str              # hotspot/scaffold/solvent/turn
    rationale: str = ""


@dataclass
class OptimizedPeptide:
    """A peptide with NCAA substitutions applied."""
    parent_name: str
    variant_name: str
    original_sequence: str
    substitutions: List[NCASubstitution]
    cyclization: str
    optimized_sequence_display: str = ""  # Human-readable with NCAA codes

    # Computed scores
    permeability_score: float = 0.0   # 0-1
    stability_score: float = 0.0      # 0-1
    affinity_retention: float = 0.0   # 0-1 (1 = full parent affinity)
    composite_score: float = 0.0
    mw_daltons: float = 0.0
    num_nme: int = 0
    num_d_aa: int = 0
    num_other_ncaa: int = 0
    drug_likeness: str = ""  # Beyond Rule of 5 assessment

    def compute_scores(self):
        """Compute optimization scores based on substitutions."""
        base_mw = len(self.original_sequence) * 110 - len(self.original_sequence) * 18
        mw_delta = sum(s.ncaa.mw_delta for s in self.substitutions)
        self.mw_daltons = round(base_mw + mw_delta, 1)

        # Permeability: base 0.2 for cyclic peptide + NCAA bonuses
        perm = 0.2
        for s in self.substitutions:
            perm += s.ncaa.permeability_bonus * 0.15  # diminishing returns
        self.permeability_score = round(min(1.0, perm), 3)

        # Stability: base 0.3 for cyclic + NCAA bonuses
        stab = 0.3
        for s in self.substitutions:
            stab += s.ncaa.stability_bonus * 0.12
        self.stability_score = round(min(1.0, stab), 3)

        # Affinity retention: start at 1.0, subtract penalties for hotspot modifications
        aff = 1.0
        for s in self.substitutions:
            if s.role == "hotspot":
                aff += s.ncaa.affinity_effect * 0.5  # hotspot changes matter more
            else:
                aff += s.ncaa.affinity_effect * 0.1  # non-hotspot changes matter less
        self.affinity_retention = round(max(0.0, min(1.0, aff)), 3)

        # Count NCAA types
        self.num_nme = sum(1 for s in self.substitutions if s.ncaa.category == "N-methyl")
        self.num_d_aa = sum(1 for s in self.substitutions if s.ncaa.category == "D-amino")
        self.num_other_ncaa = len(self.substitutions) - self.num_nme - self.num_d_aa

        # Composite score
        self.composite_score = round(
            0.35 * self.affinity_retention
            + 0.35 * self.permeability_score
            + 0.30 * self.stability_score, 3
        )

        # Drug-likeness assessment
        if self.mw_daltons <= 800:
            self.drug_likeness = "Excellent (bRo5 space)"
        elif self.mw_daltons <= 1200:
            self.drug_likeness = "Good (macrocycle space)"
        else:
            self.drug_likeness = "Challenging (biologic-like)"

    def build_display_sequence(self):
        """Build human-readable sequence with NCAA annotations."""
        seq_list = list(self.original_sequence)
        sub_map = {s.position: s.ncaa.code for s in self.substitutions}
        parts = []
        for i, aa in enumerate(seq_list):
            if i in sub_map:
                parts.append(f"[{sub_map[i]}]")
            else:
                parts.append(aa)
        self.optimized_sequence_display = "".join(parts)


# =============================================================================
# OPTIMIZATION STRATEGIES
# =============================================================================

def classify_positions(name: str, sequence: str) -> Dict[int, str]:
    """
    Classify each position in the peptide as hotspot/scaffold/solvent/turn.
    Based on parent segment information and structural knowledge.
    """
    roles = {}
    n = len(sequence)

    # Default: scaffold
    for i in range(n):
        roles[i] = "scaffold"

    if "CRD2" in name and "CRD3" not in name:
        # CRD2 designs: hotspots at positions corresponding to T62, H67, H68, Y75, Y76
        for i, aa in enumerate(sequence):
            if aa == "H" and i > 3:  # His residues are likely hotspots
                roles[i] = "hotspot"
            elif aa == "Y" and i > n // 2:  # Tyr in second half are hotspots
                roles[i] = "hotspot"
            elif aa == "T" and i < 4:  # Thr near N-term = T62 hotspot
                roles[i] = "hotspot"
            elif aa in "GPA" and i < 3:
                roles[i] = "solvent"
            elif aa == "G":
                roles[i] = "turn"
            elif aa == "P":
                roles[i] = "turn"

    elif "CRD3" in name and "CRD2" not in name:
        # CRD3 designs: hotspots at D107, N108, T110, F111
        for i, aa in enumerate(sequence):
            if aa == "D" and i < n // 2:
                roles[i] = "hotspot"
            elif aa == "N" and i < n // 2:
                roles[i] = "hotspot"
            elif aa == "T" and i > 2:
                roles[i] = "hotspot"
            elif aa == "F":
                roles[i] = "hotspot"
            elif aa in "QEL":
                roles[i] = "solvent"
            elif aa == "G":
                roles[i] = "turn"

    elif "Biloop" in name:
        # Biloop: first segment CRD2 hotspots, linker = turn, second segment CRD3 hotspots
        # Identify linker positions (G, P between segments)
        crd2_end = None
        for i, aa in enumerate(sequence):
            if i > 3 and sequence[i:i+2] in ("GP", "GG"):
                crd2_end = i
                roles[i] = "turn"
                if i + 1 < n:
                    roles[i + 1] = "turn"
                break
        # CRD2 part
        for i in range(crd2_end or n // 2):
            if sequence[i] == "H":
                roles[i] = "hotspot"
            elif sequence[i] == "T" and i == 0:
                roles[i] = "hotspot"
        # CRD3 part
        for i in range((crd2_end or n // 2) + 2, n):
            if sequence[i] in "DNTF":
                roles[i] = "hotspot"
            elif sequence[i] in "GP":
                roles[i] = "turn"

    elif "Hairpin" in name:
        # β-hairpin: Y75, Y76 in strand1; H67, H68 in strand2
        mid = n // 2
        for i, aa in enumerate(sequence):
            if aa == "Y" and i < mid:
                roles[i] = "hotspot"
            elif aa == "H" and i > mid:
                roles[i] = "hotspot"
            elif i in (mid - 1, mid):  # Turn positions
                roles[i] = "turn"
            elif aa in "GS":
                roles[i] = "solvent"

    return roles


def optimize_permeability(name: str, sequence: str, roles: Dict[int, str]) -> List[NCASubstitution]:
    """
    Apply N-methylation pattern for membrane permeability.
    Follow cyclosporine A principle: methylate ~30-50% of backbone amides,
    preferring solvent-exposed and scaffold positions.
    """
    subs = []
    n = len(sequence)
    target_nme = max(2, n // 3)  # ~33% N-methylation

    nme_map = {
        "L": "NMe-L", "G": "NMe-G", "A": "NMe-A",
        "F": "NMe-F", "T": "NMe-T", "Y": "NMe-Y",
    }

    # Priority: solvent > scaffold > turn (never methylate hotspot backbone)
    priority_order = ["solvent", "scaffold", "turn"]
    count = 0

    for priority in priority_order:
        if count >= target_nme:
            break
        for i in range(n):
            if count >= target_nme:
                break
            if roles.get(i) == priority and sequence[i] in nme_map:
                ncaa = next(nc for nc in NCAA_LIBRARY if nc.code == nme_map[sequence[i]])
                subs.append(NCASubstitution(
                    position=i,
                    original_aa=sequence[i],
                    ncaa=ncaa,
                    role=roles[i],
                    rationale=f"N-methylation at {priority} position {i} to reduce HBD count"
                ))
                count += 1

    return subs


def optimize_stability(name: str, sequence: str, roles: Dict[int, str],
                       existing_positions: set) -> List[NCASubstitution]:
    """
    Add D-amino acids and Aib at non-hotspot positions for proteolytic stability.
    """
    subs = []

    for i, aa in enumerate(sequence):
        if i in existing_positions:
            continue

        # D-Pro at turn positions
        if roles.get(i) == "turn" and aa == "P":
            ncaa = next(nc for nc in NCAA_LIBRARY if nc.code == "D-P")
            subs.append(NCASubstitution(
                position=i, original_aa=aa, ncaa=ncaa, role="turn",
                rationale="D-Pro at turn for type II' β-turn and protease resistance"
            ))
            existing_positions.add(i)

        # Aib at scaffold Ala positions
        elif roles.get(i) == "scaffold" and aa == "A" and i not in existing_positions:
            ncaa = next(nc for nc in NCAA_LIBRARY if nc.code == "Aib")
            subs.append(NCASubstitution(
                position=i, original_aa=aa, ncaa=ncaa, role="scaffold",
                rationale="Aib for backbone rigidification and protease resistance"
            ))
            existing_positions.add(i)

    return subs


def optimize_affinity(name: str, sequence: str, roles: Dict[int, str],
                      existing_positions: set) -> List[NCASubstitution]:
    """
    Enhance binding affinity at hotspot positions using conservative NCAA
    substitutions that preserve or improve interactions.
    """
    subs = []

    for i, aa in enumerate(sequence):
        if i in existing_positions:
            continue

        if roles.get(i) != "hotspot":
            continue

        # 4F-Phe at hotspot Phe positions (enhanced aromatic interaction)
        if aa == "F":
            ncaa = next(nc for nc in NCAA_LIBRARY if nc.code == "4F-Phe")
            subs.append(NCASubstitution(
                position=i, original_aa=aa, ncaa=ncaa, role="hotspot",
                rationale="4-fluoroPhe to enhance aromatic contact with TL1A"
            ))
            existing_positions.add(i)

        # Nal at Tyr hotspot (larger aromatic for enhanced π-stacking)
        elif aa == "Y" and i not in existing_positions:
            ncaa = next(nc for nc in NCAA_LIBRARY if nc.code == "Nal")
            subs.append(NCASubstitution(
                position=i, original_aa=aa, ncaa=ncaa, role="hotspot",
                rationale="Naphthylalanine for enhanced aromatic surface at hotspot"
            ))
            existing_positions.add(i)

    return subs


def generate_variants(name: str, sequence: str, cyclization: str) -> List[OptimizedPeptide]:
    """
    Generate multiple NCAA-optimized variants of a peptide design:
      Variant A: Permeability-focused (heavy N-methylation)
      Variant B: Stability-focused (D-aa + Aib)
      Variant C: Balanced (permeability + stability + affinity)
      Variant D: Maximal (all optimizations combined)
    """
    roles = classify_positions(name, sequence)
    variants = []

    # --- Variant A: Permeability ---
    subs_a = optimize_permeability(name, sequence, roles)
    var_a = OptimizedPeptide(
        parent_name=name,
        variant_name=f"{name}_permA",
        original_sequence=sequence,
        substitutions=subs_a,
        cyclization=cyclization,
    )
    var_a.compute_scores()
    var_a.build_display_sequence()
    variants.append(var_a)

    # --- Variant B: Stability ---
    subs_b = optimize_stability(name, sequence, roles, set())
    var_b = OptimizedPeptide(
        parent_name=name,
        variant_name=f"{name}_stabB",
        original_sequence=sequence,
        substitutions=subs_b,
        cyclization=cyclization,
    )
    var_b.compute_scores()
    var_b.build_display_sequence()
    variants.append(var_b)

    # --- Variant C: Balanced ---
    used_c = set()
    subs_c = optimize_permeability(name, sequence, roles)
    used_c.update(s.position for s in subs_c)
    subs_c += optimize_stability(name, sequence, roles, used_c)
    used_c.update(s.position for s in subs_c)
    subs_c += optimize_affinity(name, sequence, roles, used_c)
    var_c = OptimizedPeptide(
        parent_name=name,
        variant_name=f"{name}_balC",
        original_sequence=sequence,
        substitutions=subs_c,
        cyclization=cyclization,
    )
    var_c.compute_scores()
    var_c.build_display_sequence()
    variants.append(var_c)

    # --- Variant D: Maximal ---
    used_d = set()
    subs_d = optimize_permeability(name, sequence, roles)
    used_d.update(s.position for s in subs_d)
    subs_d += optimize_stability(name, sequence, roles, used_d)
    used_d.update(s.position for s in subs_d)
    subs_d += optimize_affinity(name, sequence, roles, used_d)
    used_d.update(s.position for s in subs_d)
    # Add extra β-amino acids at remaining scaffold positions
    for i, aa in enumerate(sequence):
        if i not in used_d and roles.get(i) == "scaffold" and aa == "F":
            ncaa = next(nc for nc in NCAA_LIBRARY if nc.code == "β3-hF")
            subs_d.append(NCASubstitution(
                position=i, original_aa=aa, ncaa=ncaa, role="scaffold",
                rationale="β3-homoPhe for maximum protease resistance at scaffold"
            ))
            used_d.add(i)
    var_d = OptimizedPeptide(
        parent_name=name,
        variant_name=f"{name}_maxD",
        original_sequence=sequence,
        substitutions=subs_d,
        cyclization=cyclization,
    )
    var_d.compute_scores()
    var_d.build_display_sequence()
    variants.append(var_d)

    return variants


# =============================================================================
# MAIN OPTIMIZATION PIPELINE
# =============================================================================

def run_optimization():
    print("=" * 70)
    print("  NCAA OPTIMIZATION — CYCLIC PEPTIDE CANDIDATES")
    print("=" * 70)

    # Load parent designs
    design_file = "output/peptide_designs.json"
    if os.path.exists(design_file):
        with open(design_file) as f:
            parent_designs = json.load(f)
    else:
        print("[!] Run 02_cyclic_peptide_design.py first to generate parent designs.")
        print("    Using built-in top candidates...")
        parent_designs = [
            {"name": "CP-Biloop-GP16", "sequence": "TCFPNHHGPDNQTFGP",
             "cyclization": "head-to-tail"},
            {"name": "CP-CRD3-9mer", "sequence": "THDNQTFPG",
             "cyclization": "head-to-tail"},
            {"name": "CP-CRD2-ext19", "sequence": "GQTCFPNHHKCDCYPGSYY",
             "cyclization": "head-to-tail"},
        ]

    # Take top candidates (by rank or first 4)
    top_parents = parent_designs[:4]

    all_variants = []
    for parent in top_parents:
        print(f"\n{'─' * 60}")
        print(f"  Parent: {parent['name']} ({parent['sequence']})")
        print(f"{'─' * 60}")

        variants = generate_variants(
            parent["name"], parent["sequence"], parent["cyclization"]
        )
        all_variants.extend(variants)

        for v in variants:
            print(f"\n  {v.variant_name}:")
            print(f"    Display:      {v.optimized_sequence_display}")
            print(f"    NCAAs:        {v.num_nme} NMe + {v.num_d_aa} D-aa "
                  f"+ {v.num_other_ncaa} other")
            print(f"    MW:           {v.mw_daltons:.0f} Da")
            print(f"    Permeability: {v.permeability_score:.3f}")
            print(f"    Stability:    {v.stability_score:.3f}")
            print(f"    Aff retain:   {v.affinity_retention:.3f}")
            print(f"    Composite:    {v.composite_score:.3f}")
            print(f"    Drug-like:    {v.drug_likeness}")
            if v.substitutions:
                print(f"    Substitutions:")
                for s in v.substitutions:
                    print(f"      pos {s.position}: {s.original_aa} → [{s.ncaa.code}] "
                          f"({s.role}) — {s.rationale}")

    # Rank all variants
    all_variants.sort(key=lambda v: v.composite_score, reverse=True)

    print(f"\n{'=' * 70}")
    print(f"  GLOBAL RANKING — ALL NCAA-OPTIMIZED VARIANTS")
    print(f"{'=' * 70}")
    print(f"\n{'Rank':<5} {'Name':<30} {'MW':<7} {'Perm':<7} {'Stab':<7} "
          f"{'Aff':<7} {'Score':<7} {'Drug-like':<20}")
    print("-" * 95)
    for i, v in enumerate(all_variants, 1):
        print(f"{i:<5} {v.variant_name:<30} {v.mw_daltons:<7.0f} "
              f"{v.permeability_score:<7.3f} {v.stability_score:<7.3f} "
              f"{v.affinity_retention:<7.3f} {v.composite_score:<7.3f} "
              f"{v.drug_likeness:<20}")

    # Save results
    os.makedirs("output", exist_ok=True)
    results = []
    for i, v in enumerate(all_variants, 1):
        entry = {
            "rank": i,
            "variant_name": v.variant_name,
            "parent_name": v.parent_name,
            "original_sequence": v.original_sequence,
            "optimized_display": v.optimized_sequence_display,
            "cyclization": v.cyclization,
            "mw_daltons": v.mw_daltons,
            "permeability_score": v.permeability_score,
            "stability_score": v.stability_score,
            "affinity_retention": v.affinity_retention,
            "composite_score": v.composite_score,
            "drug_likeness": v.drug_likeness,
            "num_nme": v.num_nme,
            "num_d_aa": v.num_d_aa,
            "num_other_ncaa": v.num_other_ncaa,
            "substitutions": [
                {
                    "position": s.position,
                    "original": s.original_aa,
                    "ncaa_code": s.ncaa.code,
                    "ncaa_name": s.ncaa.full_name,
                    "category": s.ncaa.category,
                    "role": s.role,
                    "rationale": s.rationale,
                }
                for s in v.substitutions
            ],
        }
        results.append(entry)

    with open("output/ncaa_optimized.json", "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n[✓] NCAA-optimized designs saved to output/ncaa_optimized.json")

    # Print top 3 recommendations
    print(f"\n{'=' * 70}")
    print(f"  TOP 3 RECOMMENDED CANDIDATES FOR SYNTHESIS")
    print(f"{'=' * 70}")
    for i, v in enumerate(all_variants[:3], 1):
        print(f"\n  #{i} {v.variant_name}")
        print(f"     Original:   {v.original_sequence}")
        print(f"     Optimized:  {v.optimized_sequence_display}")
        print(f"     Cyclization: {v.cyclization}")
        print(f"     MW: {v.mw_daltons:.0f} Da | Perm: {v.permeability_score:.2f} "
              f"| Stab: {v.stability_score:.2f} | Aff: {v.affinity_retention:.2f}")
        print(f"     NCAAs: {v.num_nme} NMe, {v.num_d_aa} D-aa, {v.num_other_ncaa} other")

    return all_variants


if __name__ == "__main__":
    run_optimization()
