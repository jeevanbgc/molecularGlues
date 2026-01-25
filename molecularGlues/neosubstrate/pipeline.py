"""
Neo-Substrate Generation Pipeline

This module provides the main pipeline for identifying Neo-substrates
that can form productive ternary complexes with molecular glues and
E3 ligases.

The pipeline follows a computational proteomics approach:
1. Analyze molecular glue interface
2. Generate pharmacophore model
3. Scan proteome for compatible surfaces
4. Dock top candidates
5. Score and rank results

Example:
    >>> pipeline = NeoSubstratePipeline(
    ...     e3_structure="crbn.pdb",
    ...     glue_structure="thalidomide.sdf"
    ... )
    >>> candidates = pipeline.run(
    ...     proteome_path="human_proteome/",
    ...     top_k=100
    ... )
    >>> pipeline.export_results("results.csv")
"""

import os
import json
import logging
from pathlib import Path
from dataclasses import dataclass, field, asdict
from typing import List, Dict, Optional, Tuple, Union, Callable
from datetime import datetime
import warnings

try:
    import pandas as pd
except ImportError:
    pd = None

from .interface_analysis import (
    InterfaceAnalyzer,
    InterfaceRegion,
    ContactResidue,
)
from .pharmacophore import (
    Pharmacophore,
    PharmacophoreGenerator,
    pharmacophore_to_vector,
)
from .surface_scanner import (
    SurfaceScanner,
    ProteinSurface,
    NeoSubstrateHit,
    scan_proteome_for_neo_substrates,
)
from .ternary_dock import (
    TernaryDocker,
    TernaryComplex,
    dock_neo_substrate,
    validate_ternary_geometry,
)
from .scoring import (
    NeoSubstrateScorer,
    NeoSubstrateScore,
    ScoringWeights,
)


# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@dataclass
class NeoSubstrateCandidate:
    """
    Complete representation of a Neo-substrate candidate.

    Contains all information from the pipeline including:
    - Surface patch data
    - Docking results
    - Scores and rankings
    - Validation results
    """
    protein_id: str
    structure_path: str
    patch_residues: List[Dict]
    similarity_score: float
    docking_score: Optional[float] = None
    total_score: float = 0.0
    rank: int = 0
    validation_passed: bool = False
    interface_area: float = 0.0
    pharmacophore_match: float = 0.0
    geometry_valid: bool = False
    metadata: Dict = field(default_factory=dict)

    def to_dict(self) -> Dict:
        """Convert to dictionary."""
        return asdict(self)


@dataclass
class PipelineConfig:
    """Configuration for the Neo-substrate pipeline."""
    # Surface scanning parameters
    sasa_threshold: float = 10.0
    patch_distance: float = 8.0
    min_patch_size: int = 3
    max_patch_size: int = 20
    similarity_threshold: float = 0.5

    # Docking parameters
    box_size: Tuple[float, float, float] = (30.0, 30.0, 30.0)
    exhaustiveness: int = 32
    n_poses: int = 20

    # Scoring parameters
    weights: Optional[ScoringWeights] = None

    # Pipeline parameters
    n_workers: int = 4
    dock_top_n: int = 50  # Number of candidates to dock
    validate_top_n: int = 20  # Number to validate


class NeoSubstratePipeline:
    """
    Main pipeline for Neo-substrate discovery.

    This class orchestrates the complete workflow from interface
    analysis through final ranking of candidates.

    Attributes:
        e3_structure: Path to E3 ligase structure
        glue_structure: Path to molecular glue structure
        config: Pipeline configuration
        results: List of identified candidates

    Example:
        >>> pipeline = NeoSubstratePipeline(
        ...     e3_structure="structures/crbn.pdb",
        ...     glue_structure="structures/thalidomide.sdf",
        ...     ternary_structure="structures/ternary.pdb"  # Optional
        ... )
        >>>
        >>> # Run on a set of target proteins
        >>> candidates = pipeline.run(
        ...     proteome_path="proteome/",
        ...     top_k=100
        ... )
        >>>
        >>> # Export results
        >>> pipeline.export_results("neo_substrates.csv")
    """

    def __init__(
        self,
        e3_structure: str,
        glue_structure: str,
        ternary_structure: Optional[str] = None,
        config: Optional[PipelineConfig] = None,
    ):
        """
        Initialize the pipeline.

        Args:
            e3_structure: Path to E3 ligase PDB file
            glue_structure: Path to molecular glue SDF/MOL file
            ternary_structure: Optional path to existing ternary complex
            config: Pipeline configuration
        """
        self.e3_structure = Path(e3_structure)
        self.glue_structure = Path(glue_structure)
        self.ternary_structure = Path(ternary_structure) if ternary_structure else None
        self.config = config or PipelineConfig()

        # Initialize components
        self._interface_analyzer = None
        self._pharmacophore_gen = PharmacophoreGenerator()
        self._surface_scanner = None
        self._docker = None
        self._scorer = None

        # Results storage
        self.interface: Optional[InterfaceRegion] = None
        self.reference_pharmacophore: Optional[Pharmacophore] = None
        self.hits: List[NeoSubstrateHit] = []
        self.candidates: List[NeoSubstrateCandidate] = []

        # Callbacks
        self._progress_callback: Optional[Callable] = None

    def set_progress_callback(self, callback: Callable[[str, int, int], None]):
        """
        Set progress callback for pipeline steps.

        Args:
            callback: Function(step_name, current, total)
        """
        self._progress_callback = callback

    def _report_progress(self, step: str, current: int, total: int):
        """Report progress to callback if set."""
        if self._progress_callback:
            self._progress_callback(step, current, total)
        logger.info(f"{step}: {current}/{total}")

    def analyze_interface(self) -> InterfaceRegion:
        """
        Analyze the molecular glue binding interface.

        Returns:
            InterfaceRegion with contact analysis
        """
        logger.info("Analyzing molecular glue interface...")

        if self.ternary_structure:
            # Analyze existing ternary complex
            self._interface_analyzer = InterfaceAnalyzer(
                structure_path=str(self.ternary_structure),
                glue_selection="resname LIG or resname UNK",
            )
        else:
            # Analyze E3-glue binary complex
            self._interface_analyzer = InterfaceAnalyzer(
                structure_path=str(self.e3_structure),
                glue_selection="resname LIG or resname UNK",
            )

        self.interface = self._interface_analyzer.analyze()

        logger.info(
            f"Found {len(self.interface.e3_contacts)} E3 contacts, "
            f"{len(self.interface.substrate_contacts)} substrate contacts"
        )

        return self.interface

    def generate_pharmacophore(self) -> Pharmacophore:
        """
        Generate reference pharmacophore from interface.

        Returns:
            Pharmacophore model for Neo-substrate matching
        """
        logger.info("Generating reference pharmacophore...")

        if self.interface is None:
            self.analyze_interface()

        # Convert contacts to format expected by pharmacophore generator
        contacts = []

        # Use substrate contacts if available (from ternary)
        # Otherwise use E3 contacts as proxy
        source_contacts = (
            self.interface.substrate_contacts
            if self.interface.substrate_contacts
            else self.interface.e3_contacts
        )

        for contact in source_contacts:
            contacts.append({
                'resid': contact.resid,
                'resname': contact.resname,
                'chain': contact.chain,
                'min_distance': contact.min_distance,
            })

        self.reference_pharmacophore = self._pharmacophore_gen.from_interface_contacts(
            contacts
        )
        self.reference_pharmacophore.name = "glue_interface"

        logger.info(
            f"Generated pharmacophore with {len(self.reference_pharmacophore.features)} features"
        )

        return self.reference_pharmacophore

    def scan_proteome(
        self,
        proteome_path: str,
        file_pattern: str = "*.pdb",
        similarity_threshold: Optional[float] = None,
    ) -> List[NeoSubstrateHit]:
        """
        Scan proteome for compatible Neo-substrates.

        Args:
            proteome_path: Directory containing protein structures
            file_pattern: Glob pattern for structure files
            similarity_threshold: Minimum similarity for hits

        Returns:
            List of NeoSubstrateHit objects
        """
        logger.info(f"Scanning proteome at {proteome_path}...")

        if self.reference_pharmacophore is None:
            self.generate_pharmacophore()

        threshold = similarity_threshold or self.config.similarity_threshold

        # Find all structure files
        proteome_dir = Path(proteome_path)
        structure_files = list(proteome_dir.glob(file_pattern))

        if not structure_files:
            logger.warning(f"No structure files found matching {file_pattern}")
            return []

        logger.info(f"Found {len(structure_files)} structures to scan")

        # Scan proteome
        def progress_cb(current, total):
            self._report_progress("Scanning", current, total)

        self.hits = scan_proteome_for_neo_substrates(
            glue_pharmacophore=self.reference_pharmacophore,
            proteome_structures=[str(f) for f in structure_files],
            similarity_threshold=threshold,
            n_workers=self.config.n_workers,
            progress_callback=progress_cb,
        )

        logger.info(f"Found {len(self.hits)} potential Neo-substrate hits")

        return self.hits

    def dock_candidates(
        self,
        top_n: Optional[int] = None,
    ) -> Dict[str, TernaryComplex]:
        """
        Dock top candidates to form ternary complexes.

        Args:
            top_n: Number of top candidates to dock

        Returns:
            Dictionary of protein_id -> TernaryComplex
        """
        n = top_n or self.config.dock_top_n

        if not self.hits:
            logger.warning("No hits to dock. Run scan_proteome first.")
            return {}

        # Take top N hits
        to_dock = self.hits[:n]
        logger.info(f"Docking top {len(to_dock)} candidates...")

        self._docker = TernaryDocker(
            e3_pdb=str(self.e3_structure),
            glue_sdf=str(self.glue_structure),
            box_size=self.config.box_size,
            exhaustiveness=self.config.exhaustiveness,
            n_poses=self.config.n_poses,
        )

        docking_results = {}

        for i, hit in enumerate(to_dock):
            self._report_progress("Docking", i + 1, len(to_dock))

            try:
                complex_result = self._docker.dock_substrate(hit.structure_path)
                docking_results[hit.protein_id] = complex_result
            except Exception as e:
                logger.warning(f"Docking failed for {hit.protein_id}: {e}")

        logger.info(f"Successfully docked {len(docking_results)} candidates")

        return docking_results

    def score_and_rank(
        self,
        docking_results: Optional[Dict[str, TernaryComplex]] = None,
    ) -> List[NeoSubstrateCandidate]:
        """
        Score and rank all candidates.

        Args:
            docking_results: Optional docking results

        Returns:
            List of ranked NeoSubstrateCandidate objects
        """
        logger.info("Scoring and ranking candidates...")

        self._scorer = NeoSubstrateScorer(
            reference_pharmacophore=self.reference_pharmacophore,
            weights=self.config.weights,
        )

        # Score candidates
        scores = self._scorer.score_candidates(self.hits, docking_results)

        # Rank
        ranked_scores = self._scorer.rank(scores)

        # Convert to NeoSubstrateCandidate
        self.candidates = []

        for score in ranked_scores:
            # Find corresponding hit
            hit = next(
                (h for h in self.hits if h.protein_id == score.candidate_id),
                None
            )

            if hit is None:
                continue

            # Get docking score if available
            docking_score = None
            if docking_results and score.candidate_id in docking_results:
                docking_score = docking_results[score.candidate_id].total_score

            candidate = NeoSubstrateCandidate(
                protein_id=score.candidate_id,
                structure_path=hit.structure_path,
                patch_residues=[
                    {
                        'resid': r.resid,
                        'resname': r.resname,
                        'chain': r.chain,
                    }
                    for r in hit.patch.residues
                ],
                similarity_score=hit.similarity_score,
                docking_score=docking_score,
                total_score=score.total_score,
                rank=score.rank,
                interface_area=hit.patch.area,
                pharmacophore_match=score.pharmacophore_score,
                metadata={
                    'scores': score.to_dict(),
                }
            )

            self.candidates.append(candidate)

        logger.info(f"Ranked {len(self.candidates)} candidates")

        return self.candidates

    def validate_top_candidates(
        self,
        top_n: Optional[int] = None,
    ) -> List[NeoSubstrateCandidate]:
        """
        Validate geometry of top candidates.

        Args:
            top_n: Number of top candidates to validate

        Returns:
            List of validated candidates
        """
        n = top_n or self.config.validate_top_n

        to_validate = self.candidates[:n]
        logger.info(f"Validating top {len(to_validate)} candidates...")

        for i, candidate in enumerate(to_validate):
            self._report_progress("Validating", i + 1, len(to_validate))

            try:
                validation = validate_ternary_geometry(
                    e3_pdb=str(self.e3_structure),
                    substrate_pdb=candidate.structure_path,
                    glue_sdf=str(self.glue_structure),
                )

                candidate.validation_passed = validation['is_valid']
                candidate.geometry_valid = validation['is_valid']
                candidate.metadata['validation'] = validation

            except Exception as e:
                logger.warning(f"Validation failed for {candidate.protein_id}: {e}")
                candidate.validation_passed = False

        validated = [c for c in to_validate if c.validation_passed]
        logger.info(f"{len(validated)} candidates passed validation")

        return validated

    def run(
        self,
        proteome_path: str,
        file_pattern: str = "*.pdb",
        top_k: int = 100,
        dock: bool = True,
        validate: bool = True,
    ) -> List[NeoSubstrateCandidate]:
        """
        Run the complete Neo-substrate discovery pipeline.

        Args:
            proteome_path: Directory containing protein structures
            file_pattern: Glob pattern for structure files
            top_k: Number of top candidates to return
            dock: Whether to perform docking
            validate: Whether to validate geometry

        Returns:
            List of top NeoSubstrateCandidate objects
        """
        logger.info("=" * 50)
        logger.info("Starting Neo-Substrate Discovery Pipeline")
        logger.info("=" * 50)

        start_time = datetime.now()

        # Step 1: Analyze interface
        self.analyze_interface()

        # Step 2: Generate pharmacophore
        self.generate_pharmacophore()

        # Step 3: Scan proteome
        self.scan_proteome(proteome_path, file_pattern)

        # Step 4: Dock candidates (optional)
        docking_results = None
        if dock and self.hits:
            docking_results = self.dock_candidates()

        # Step 5: Score and rank
        self.score_and_rank(docking_results)

        # Step 6: Validate (optional)
        if validate and self.candidates:
            self.validate_top_candidates()

        elapsed = datetime.now() - start_time
        logger.info(f"Pipeline completed in {elapsed}")
        logger.info(f"Found {len(self.candidates)} candidates")

        return self.candidates[:top_k]

    def export_results(
        self,
        output_path: str,
        format: str = "csv",
    ):
        """
        Export results to file.

        Args:
            output_path: Output file path
            format: Output format ("csv", "json", "tsv")
        """
        if not self.candidates:
            logger.warning("No candidates to export")
            return

        output_path = Path(output_path)

        if format == "json":
            data = [c.to_dict() for c in self.candidates]
            with open(output_path, 'w') as f:
                json.dump(data, f, indent=2, default=str)

        elif format in ("csv", "tsv"):
            if pd is None:
                # Fallback without pandas
                delimiter = "\t" if format == "tsv" else ","
                with open(output_path, 'w') as f:
                    # Header
                    headers = [
                        "rank", "protein_id", "total_score", "similarity_score",
                        "docking_score", "interface_area", "validation_passed"
                    ]
                    f.write(delimiter.join(headers) + "\n")

                    # Data
                    for c in self.candidates:
                        row = [
                            str(c.rank),
                            c.protein_id,
                            f"{c.total_score:.4f}",
                            f"{c.similarity_score:.4f}",
                            f"{c.docking_score:.2f}" if c.docking_score else "N/A",
                            f"{c.interface_area:.1f}",
                            str(c.validation_passed),
                        ]
                        f.write(delimiter.join(row) + "\n")
            else:
                # Use pandas
                data = [c.to_dict() for c in self.candidates]
                df = pd.DataFrame(data)
                if format == "csv":
                    df.to_csv(output_path, index=False)
                else:
                    df.to_csv(output_path, index=False, sep="\t")

        logger.info(f"Exported {len(self.candidates)} candidates to {output_path}")

    def get_summary(self) -> Dict:
        """
        Get pipeline summary statistics.

        Returns:
            Dictionary with summary statistics
        """
        return {
            'e3_structure': str(self.e3_structure),
            'glue_structure': str(self.glue_structure),
            'interface_contacts': (
                len(self.interface.e3_contacts) if self.interface else 0
            ),
            'pharmacophore_features': (
                len(self.reference_pharmacophore.features)
                if self.reference_pharmacophore else 0
            ),
            'total_hits': len(self.hits),
            'total_candidates': len(self.candidates),
            'validated_candidates': sum(
                1 for c in self.candidates if c.validation_passed
            ),
            'top_score': (
                self.candidates[0].total_score if self.candidates else 0
            ),
        }


def run_quick_scan(
    e3_structure: str,
    glue_structure: str,
    target_structures: List[str],
    similarity_threshold: float = 0.5,
) -> List[Dict]:
    """
    Run a quick scan without docking.

    Useful for initial filtering before more expensive calculations.

    Args:
        e3_structure: Path to E3 ligase structure
        glue_structure: Path to molecular glue structure
        target_structures: List of paths to target structures
        similarity_threshold: Minimum similarity

    Returns:
        List of hit dictionaries
    """
    pipeline = NeoSubstratePipeline(
        e3_structure=e3_structure,
        glue_structure=glue_structure,
    )

    # Generate reference pharmacophore
    pipeline.analyze_interface()
    pipeline.generate_pharmacophore()

    # Scan targets
    hits = scan_proteome_for_neo_substrates(
        glue_pharmacophore=pipeline.reference_pharmacophore,
        proteome_structures=target_structures,
        similarity_threshold=similarity_threshold,
        n_workers=1,
    )

    return [
        {
            'protein_id': h.protein_id,
            'similarity_score': h.similarity_score,
            'structure_path': h.structure_path,
            'patch_area': h.patch.area,
            'n_residues': len(h.patch.residues),
        }
        for h in hits
    ]
