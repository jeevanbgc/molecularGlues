"""
Neo-Substrate Generation Module for Molecular Glue Ternary Complexes

This module provides computational tools for identifying and validating
Neo-substrates that can form productive ternary complexes with molecular
glues and E3 ligases.

Workflow:
    1. Interface Analysis - Characterize molecular glue binding interface
    2. Pharmacophore Generation - Extract key interaction features
    3. Surface Scanning - Screen proteome for compatible surfaces
    4. Ternary Docking - Validate candidate Neo-substrates
    5. Scoring & Ranking - Prioritize candidates for experimental validation

Example:
    >>> from neosubstrate import NeoSubstratePipeline
    >>> pipeline = NeoSubstratePipeline(
    ...     e3_structure="crbn.pdb",
    ...     glue_structure="thalidomide.sdf"
    ... )
    >>> candidates = pipeline.run(proteome_db="human_proteome.fasta")
"""

__version__ = "0.1.0"
__author__ = "Molecular Glues Project"

from .interface_analysis import (
    InterfaceAnalyzer,
    identify_glue_interface,
    calculate_interface_area,
    get_contact_residues,
)

from .pharmacophore import (
    PharmacophoreGenerator,
    pharmacophore_to_vector,
    calculate_pharmacophore_similarity,
)

from .surface_scanner import (
    SurfaceScanner,
    ProteinSurface,
    scan_proteome_for_neo_substrates,
)

from .ternary_dock import (
    TernaryDocker,
    dock_neo_substrate,
    validate_ternary_geometry,
)

from .scoring import (
    NeoSubstrateScorer,
    calculate_ternary_score,
    rank_candidates,
)

from .pipeline import (
    NeoSubstratePipeline,
    NeoSubstrateCandidate,
)

__all__ = [
    # Interface Analysis
    "InterfaceAnalyzer",
    "identify_glue_interface",
    "calculate_interface_area",
    "get_contact_residues",
    # Pharmacophore
    "PharmacophoreGenerator",
    "pharmacophore_to_vector",
    "calculate_pharmacophore_similarity",
    # Surface Scanner
    "SurfaceScanner",
    "ProteinSurface",
    "scan_proteome_for_neo_substrates",
    # Ternary Docking
    "TernaryDocker",
    "dock_neo_substrate",
    "validate_ternary_geometry",
    # Scoring
    "NeoSubstrateScorer",
    "calculate_ternary_score",
    "rank_candidates",
    # Pipeline
    "NeoSubstratePipeline",
    "NeoSubstrateCandidate",
]
