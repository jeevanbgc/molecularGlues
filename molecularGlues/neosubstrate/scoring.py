"""
Scoring Module for Neo-Substrate Ranking

This module provides comprehensive scoring functions for ranking
Neo-substrate candidates based on multiple criteria including:

- Pharmacophore similarity
- Docking scores
- Interface quality
- Druggability metrics
- Structural features

The scoring combines physics-based and knowledge-based terms
to prioritize candidates for experimental validation.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional, Union
from pathlib import Path
import warnings
from enum import Enum

try:
    from sklearn.preprocessing import StandardScaler, MinMaxScaler
    from sklearn.ensemble import RandomForestRegressor
except ImportError:
    StandardScaler = None
    RandomForestRegressor = None

try:
    import tensorflow as tf
    from tensorflow import keras
except ImportError:
    tf = None
    keras = None

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, AllChem, rdMolDescriptors
    from rdkit.Chem import MACCSkeys
except ImportError:
    Chem = None

from .pharmacophore import (
    Pharmacophore,
    pharmacophore_to_vector,
    calculate_pharmacophore_similarity,
)
from .surface_scanner import SurfacePatch, NeoSubstrateHit
from .ternary_dock import TernaryComplex, DockingPose


class ScoreComponent(Enum):
    """Components of the Neo-substrate score."""
    PHARMACOPHORE_SIMILARITY = "pharmacophore_similarity"
    DOCKING_SCORE = "docking_score"
    INTERFACE_QUALITY = "interface_quality"
    DRUGGABILITY = "druggability"
    STRUCTURAL_FEATURES = "structural_features"
    DEGRADATION_POTENTIAL = "degradation_potential"


@dataclass
class ScoringWeights:
    """Weights for different scoring components."""
    pharmacophore: float = 0.25
    docking: float = 0.25
    interface: float = 0.20
    druggability: float = 0.15
    structural: float = 0.10
    degradation: float = 0.05

    def normalize(self):
        """Normalize weights to sum to 1.0."""
        total = (
            self.pharmacophore + self.docking + self.interface +
            self.druggability + self.structural + self.degradation
        )
        if total > 0:
            self.pharmacophore /= total
            self.docking /= total
            self.interface /= total
            self.druggability /= total
            self.structural /= total
            self.degradation /= total


@dataclass
class NeoSubstrateScore:
    """
    Comprehensive score for a Neo-substrate candidate.

    Contains individual component scores and the final weighted score.
    """
    candidate_id: str
    total_score: float = 0.0
    pharmacophore_score: float = 0.0
    docking_score: float = 0.0
    interface_score: float = 0.0
    druggability_score: float = 0.0
    structural_score: float = 0.0
    degradation_score: float = 0.0
    rank: int = 0
    confidence: float = 0.0
    details: Dict = field(default_factory=dict)

    def to_dict(self) -> Dict:
        """Convert to dictionary."""
        return {
            'candidate_id': self.candidate_id,
            'total_score': self.total_score,
            'rank': self.rank,
            'confidence': self.confidence,
            'components': {
                'pharmacophore': self.pharmacophore_score,
                'docking': self.docking_score,
                'interface': self.interface_score,
                'druggability': self.druggability_score,
                'structural': self.structural_score,
                'degradation': self.degradation_score,
            },
            'details': self.details
        }


class NeoSubstrateScorer:
    """
    Comprehensive scorer for Neo-substrate candidates.

    This class combines multiple scoring functions to rank
    Neo-substrate candidates for molecular glue-mediated
    degradation.

    Example:
        >>> scorer = NeoSubstrateScorer(reference_pharmacophore=ref_pharm)
        >>> scores = scorer.score_candidates(candidates)
        >>> ranked = scorer.rank(scores)
    """

    def __init__(
        self,
        reference_pharmacophore: Optional[Pharmacophore] = None,
        weights: Optional[ScoringWeights] = None,
        normalize_scores: bool = True,
    ):
        """
        Initialize the scorer.

        Args:
            reference_pharmacophore: Reference pharmacophore for similarity
            weights: Scoring weights for different components
            normalize_scores: Whether to normalize component scores
        """
        self.reference_pharmacophore = reference_pharmacophore
        self.weights = weights or ScoringWeights()
        self.weights.normalize()
        self.normalize_scores = normalize_scores

        self._scaler = StandardScaler() if StandardScaler else None

    def score_pharmacophore_similarity(
        self,
        candidate_pharmacophore: Pharmacophore,
        method: str = "cosine"
    ) -> float:
        """
        Score pharmacophore similarity to reference.

        Args:
            candidate_pharmacophore: Candidate's pharmacophore
            method: Similarity method

        Returns:
            Similarity score (0-1)
        """
        if self.reference_pharmacophore is None:
            return 0.5  # Neutral score

        return calculate_pharmacophore_similarity(
            self.reference_pharmacophore,
            candidate_pharmacophore,
            method=method
        )

    def score_docking(
        self,
        docking_result: Union[TernaryComplex, float],
        best_known: float = -10.0,
        worst_acceptable: float = -2.0
    ) -> float:
        """
        Convert docking score to normalized score.

        Args:
            docking_result: TernaryComplex or raw docking score
            best_known: Best known docking score (kcal/mol)
            worst_acceptable: Worst acceptable score

        Returns:
            Normalized score (0-1)
        """
        if isinstance(docking_result, TernaryComplex):
            score = docking_result.total_score
        else:
            score = float(docking_result)

        # Normalize: better (more negative) scores -> higher values
        if score <= best_known:
            return 1.0
        elif score >= worst_acceptable:
            return 0.0
        else:
            return (worst_acceptable - score) / (worst_acceptable - best_known)

    def score_interface_quality(
        self,
        patch: SurfacePatch,
        min_area: float = 200.0,
        optimal_area: float = 800.0
    ) -> float:
        """
        Score the quality of the binding interface.

        Args:
            patch: Surface patch representing interface
            min_area: Minimum acceptable area (Å²)
            optimal_area: Optimal interface area

        Returns:
            Interface quality score (0-1)
        """
        area = patch.area

        if area < min_area:
            return area / min_area * 0.5

        if area <= optimal_area:
            return 0.5 + 0.5 * (area - min_area) / (optimal_area - min_area)

        # Penalty for very large interfaces (may be non-specific)
        return max(0.5, 1.0 - 0.1 * (area - optimal_area) / optimal_area)

    def score_druggability(
        self,
        structure_path: Optional[str] = None,
        patch: Optional[SurfacePatch] = None
    ) -> float:
        """
        Score the druggability of the binding site.

        Based on pocket properties like:
        - Hydrophobicity
        - Enclosure
        - Feature diversity

        Args:
            structure_path: Path to protein structure
            patch: Surface patch

        Returns:
            Druggability score (0-1)
        """
        if patch is None:
            return 0.5  # Neutral

        # Score based on feature diversity
        features = patch.get_feature_counts()

        # Prefer diverse feature composition
        non_zero = sum(1 for v in features.values() if v > 0)
        diversity = non_zero / len(features)

        # Prefer balanced hydrophobic/polar
        total_features = sum(features.values())
        if total_features == 0:
            return 0.3

        from .pharmacophore import FeatureType
        hydrophobic = features.get(FeatureType.HYDROPHOBIC, 0)
        polar = (
            features.get(FeatureType.HBOND_DONOR, 0) +
            features.get(FeatureType.HBOND_ACCEPTOR, 0)
        )

        ratio = hydrophobic / total_features
        balance = 1.0 - abs(ratio - 0.5) * 2  # Peak at 50% hydrophobic

        return 0.5 * diversity + 0.5 * balance

    def score_structural_features(
        self,
        protein_id: str,
        structure_path: Optional[str] = None
    ) -> float:
        """
        Score structural features relevant to degradation.

        Considers:
        - Lysine accessibility for ubiquitination
        - Surface exposure
        - Flexibility indicators

        Args:
            protein_id: Protein identifier
            structure_path: Path to structure file

        Returns:
            Structural features score (0-1)
        """
        # Placeholder - would need actual structure analysis
        # In practice, would check:
        # 1. Surface lysine count and accessibility
        # 2. Distance from binding site to nearest lysine
        # 3. Loop regions that might facilitate ubiquitin transfer

        return 0.5  # Neutral default

    def score_degradation_potential(
        self,
        candidate_id: str,
        known_substrates: Optional[List[str]] = None
    ) -> float:
        """
        Score based on known degradation data.

        Args:
            candidate_id: Candidate protein identifier
            known_substrates: List of known degradable substrates

        Returns:
            Degradation potential score (0-1)
        """
        if known_substrates is None:
            return 0.5  # No prior information

        # Check if candidate is similar to known substrates
        # Placeholder - would need sequence/structure comparison

        return 0.5

    def score_candidate(
        self,
        hit: NeoSubstrateHit,
        docking_result: Optional[TernaryComplex] = None
    ) -> NeoSubstrateScore:
        """
        Compute comprehensive score for a candidate.

        Args:
            hit: NeoSubstrateHit from surface scanning
            docking_result: Optional docking result

        Returns:
            NeoSubstrateScore with all components
        """
        score = NeoSubstrateScore(candidate_id=hit.protein_id)

        # Pharmacophore similarity
        patch_pharmacophore = hit.patch.to_pharmacophore()
        score.pharmacophore_score = self.score_pharmacophore_similarity(
            patch_pharmacophore
        )

        # Docking score
        if docking_result is not None:
            score.docking_score = self.score_docking(docking_result)
        else:
            score.docking_score = 0.5  # Neutral if not docked

        # Interface quality
        score.interface_score = self.score_interface_quality(hit.patch)

        # Druggability
        score.druggability_score = self.score_druggability(
            structure_path=hit.structure_path,
            patch=hit.patch
        )

        # Structural features
        score.structural_score = self.score_structural_features(
            hit.protein_id,
            hit.structure_path
        )

        # Degradation potential
        score.degradation_score = self.score_degradation_potential(
            hit.protein_id
        )

        # Calculate weighted total
        score.total_score = (
            self.weights.pharmacophore * score.pharmacophore_score +
            self.weights.docking * score.docking_score +
            self.weights.interface * score.interface_score +
            self.weights.druggability * score.druggability_score +
            self.weights.structural * score.structural_score +
            self.weights.degradation * score.degradation_score
        )

        # Estimate confidence based on data availability
        data_points = sum([
            1 if score.pharmacophore_score != 0.5 else 0,
            1 if score.docking_score != 0.5 else 0,
            1 if score.interface_score != 0.5 else 0,
            1 if score.druggability_score != 0.5 else 0,
        ])
        score.confidence = data_points / 4.0

        return score

    def score_candidates(
        self,
        hits: List[NeoSubstrateHit],
        docking_results: Optional[Dict[str, TernaryComplex]] = None
    ) -> List[NeoSubstrateScore]:
        """
        Score multiple candidates.

        Args:
            hits: List of NeoSubstrateHit objects
            docking_results: Optional dict of protein_id -> TernaryComplex

        Returns:
            List of NeoSubstrateScore objects
        """
        scores = []
        docking_results = docking_results or {}

        for hit in hits:
            docking = docking_results.get(hit.protein_id)
            score = self.score_candidate(hit, docking)
            scores.append(score)

        return scores

    def rank(self, scores: List[NeoSubstrateScore]) -> List[NeoSubstrateScore]:
        """
        Rank scored candidates.

        Args:
            scores: List of NeoSubstrateScore objects

        Returns:
            Sorted list with ranks assigned
        """
        # Sort by total score (descending)
        sorted_scores = sorted(
            scores,
            key=lambda x: x.total_score,
            reverse=True
        )

        # Assign ranks
        for i, score in enumerate(sorted_scores):
            score.rank = i + 1

        return sorted_scores


def calculate_ternary_score(
    pharmacophore_sim: float,
    docking_score: float,
    interface_area: float,
    weights: Optional[ScoringWeights] = None
) -> float:
    """
    Calculate a quick ternary complex score.

    Args:
        pharmacophore_sim: Pharmacophore similarity (0-1)
        docking_score: Docking score (kcal/mol)
        interface_area: Interface area (Å²)
        weights: Optional scoring weights

    Returns:
        Combined score (0-1)
    """
    weights = weights or ScoringWeights()

    # Normalize docking score
    normalized_docking = max(0, min(1, (-docking_score + 2) / 12))

    # Normalize interface area
    normalized_interface = min(1, interface_area / 800)

    return (
        weights.pharmacophore * pharmacophore_sim +
        weights.docking * normalized_docking +
        weights.interface * normalized_interface
    ) / (weights.pharmacophore + weights.docking + weights.interface)


def rank_candidates(
    candidates: List[Dict],
    reference_pharmacophore: Optional[Pharmacophore] = None
) -> List[Dict]:
    """
    Quick ranking function for candidate dictionaries.

    Args:
        candidates: List of candidate dictionaries
        reference_pharmacophore: Reference pharmacophore

    Returns:
        Sorted list with ranks
    """
    # Score each candidate
    for cand in candidates:
        score = 0.0

        if 'similarity_score' in cand:
            score += 0.4 * cand['similarity_score']

        if 'docking_score' in cand:
            # Normalize docking score
            dock = cand['docking_score']
            normalized = max(0, min(1, (-dock + 2) / 12))
            score += 0.4 * normalized

        if 'interface_area' in cand:
            area = cand['interface_area']
            normalized = min(1, area / 800)
            score += 0.2 * normalized

        cand['total_score'] = score

    # Sort by score
    sorted_candidates = sorted(
        candidates,
        key=lambda x: x.get('total_score', 0),
        reverse=True
    )

    # Assign ranks
    for i, cand in enumerate(sorted_candidates):
        cand['rank'] = i + 1

    return sorted_candidates


class MLScorer:
    """
    Machine learning-based scorer for Neo-substrates.

    Uses trained models to predict degradation efficiency
    based on structural and chemical features.
    """

    def __init__(self, model_path: Optional[str] = None):
        """
        Initialize the ML scorer.

        Args:
            model_path: Path to pre-trained model
        """
        self.model = None
        self.scaler = StandardScaler() if StandardScaler else None

        if model_path:
            self.load_model(model_path)

    def load_model(self, model_path: str):
        """Load a pre-trained model."""
        if keras is not None:
            try:
                self.model = keras.models.load_model(model_path)
            except Exception:
                pass

    def featurize(
        self,
        hit: NeoSubstrateHit,
        docking_result: Optional[TernaryComplex] = None
    ) -> np.ndarray:
        """
        Convert a hit to feature vector.

        Args:
            hit: NeoSubstrateHit object
            docking_result: Optional docking result

        Returns:
            Feature vector
        """
        features = []

        # Pharmacophore features
        pharma = hit.patch.to_pharmacophore()
        pharma_vec = pharmacophore_to_vector(pharma)
        features.extend(pharma_vec.tolist())

        # Interface features
        features.append(hit.similarity_score)
        features.append(hit.patch.area)
        features.append(len(hit.patch.residues))

        # Docking features
        if docking_result:
            features.append(docking_result.total_score)
        else:
            features.append(0.0)

        return np.array(features, dtype=np.float32)

    def predict(
        self,
        hits: List[NeoSubstrateHit],
        docking_results: Optional[Dict[str, TernaryComplex]] = None
    ) -> List[float]:
        """
        Predict degradation scores for candidates.

        Args:
            hits: List of NeoSubstrateHit objects
            docking_results: Optional docking results

        Returns:
            List of predicted scores
        """
        if self.model is None:
            warnings.warn("No model loaded, returning similarity scores")
            return [hit.similarity_score for hit in hits]

        docking_results = docking_results or {}

        # Featurize
        features = np.array([
            self.featurize(hit, docking_results.get(hit.protein_id))
            for hit in hits
        ])

        # Scale
        if self.scaler:
            features = self.scaler.fit_transform(features)

        # Predict
        predictions = self.model.predict(features)

        return predictions.flatten().tolist()

    def build_model(self, input_dim: int = 32) -> 'keras.Model':
        """
        Build a neural network model for scoring.

        Args:
            input_dim: Input feature dimension

        Returns:
            Keras model
        """
        if keras is None:
            raise ImportError("TensorFlow/Keras required for ML scoring")

        model = keras.Sequential([
            keras.layers.Dense(64, activation='relu', input_shape=(input_dim,)),
            keras.layers.Dropout(0.3),
            keras.layers.Dense(32, activation='relu'),
            keras.layers.Dropout(0.2),
            keras.layers.Dense(16, activation='relu'),
            keras.layers.Dense(1, activation='sigmoid')
        ])

        model.compile(
            optimizer='adam',
            loss='binary_crossentropy',
            metrics=['accuracy']
        )

        self.model = model
        return model
