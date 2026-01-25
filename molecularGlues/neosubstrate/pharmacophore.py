"""
Pharmacophore Generation Module for Neo-Substrate Discovery

This module extracts pharmacophore features from molecular glue interfaces
and uses them to identify compatible Neo-substrate surfaces. The approach
is analogous to computational proteomics methods for identifying protein-protein
interaction surfaces.

Key Features:
    - Extract 3D pharmacophore from interface contacts
    - Generate feature vectors for similarity comparison
    - Support for multiple pharmacophore types (hydrophobic, HBD, HBA, etc.)
    - Pharmacophore-based virtual screening
"""

import numpy as np
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional, Set, Union
from pathlib import Path
import warnings
from enum import Enum

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors
    from rdkit.Chem.Pharm3D import Pharmacophore as RDPharm
    from rdkit.Chem.Pharm3D import EmbedLib
    from rdkit.Chem import ChemicalFeatures
    from rdkit import RDConfig
    import os
    FDEF_FILE = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
except ImportError:
    Chem = None
    FDEF_FILE = None
    warnings.warn("RDKit not installed. Pharmacophore features will be limited.")

try:
    from sklearn.metrics.pairwise import cosine_similarity, euclidean_distances
    from sklearn.preprocessing import StandardScaler
except ImportError:
    cosine_similarity = None
    warnings.warn("scikit-learn not installed. Similarity calculations will be limited.")


class FeatureType(Enum):
    """Pharmacophore feature types based on interaction chemistry."""
    HYDROPHOBIC = "hydrophobic"
    AROMATIC = "aromatic"
    HBOND_DONOR = "hbond_donor"
    HBOND_ACCEPTOR = "hbond_acceptor"
    POSITIVE = "positive_ionizable"
    NEGATIVE = "negative_ionizable"
    HALOGEN = "halogen"


# Mapping of residue names to pharmacophore features
RESIDUE_FEATURES = {
    # Hydrophobic residues
    'ALA': [FeatureType.HYDROPHOBIC],
    'VAL': [FeatureType.HYDROPHOBIC],
    'LEU': [FeatureType.HYDROPHOBIC],
    'ILE': [FeatureType.HYDROPHOBIC],
    'MET': [FeatureType.HYDROPHOBIC],
    'PRO': [FeatureType.HYDROPHOBIC],
    # Aromatic residues
    'PHE': [FeatureType.HYDROPHOBIC, FeatureType.AROMATIC],
    'TYR': [FeatureType.HYDROPHOBIC, FeatureType.AROMATIC, FeatureType.HBOND_DONOR, FeatureType.HBOND_ACCEPTOR],
    'TRP': [FeatureType.HYDROPHOBIC, FeatureType.AROMATIC, FeatureType.HBOND_DONOR],
    # Polar residues
    'SER': [FeatureType.HBOND_DONOR, FeatureType.HBOND_ACCEPTOR],
    'THR': [FeatureType.HBOND_DONOR, FeatureType.HBOND_ACCEPTOR],
    'CYS': [FeatureType.HBOND_DONOR],
    'ASN': [FeatureType.HBOND_DONOR, FeatureType.HBOND_ACCEPTOR],
    'GLN': [FeatureType.HBOND_DONOR, FeatureType.HBOND_ACCEPTOR],
    # Charged residues
    'LYS': [FeatureType.POSITIVE, FeatureType.HBOND_DONOR],
    'ARG': [FeatureType.POSITIVE, FeatureType.HBOND_DONOR],
    'HIS': [FeatureType.POSITIVE, FeatureType.HBOND_DONOR, FeatureType.HBOND_ACCEPTOR, FeatureType.AROMATIC],
    'ASP': [FeatureType.NEGATIVE, FeatureType.HBOND_ACCEPTOR],
    'GLU': [FeatureType.NEGATIVE, FeatureType.HBOND_ACCEPTOR],
    # Special
    'GLY': [],  # Too small for significant interactions
}


@dataclass
class PharmacophoreFeature:
    """
    Represents a single pharmacophore feature point.

    Attributes:
        feature_type: Type of pharmacophore feature
        position: 3D coordinates (x, y, z) in Angstroms
        radius: Tolerance radius for feature matching
        source: Origin of feature (residue ID or atom name)
        weight: Importance weight for scoring
    """
    feature_type: FeatureType
    position: Tuple[float, float, float]
    radius: float = 1.5
    source: str = ""
    weight: float = 1.0

    def distance_to(self, other: 'PharmacophoreFeature') -> float:
        """Calculate Euclidean distance to another feature."""
        return np.sqrt(sum(
            (a - b) ** 2
            for a, b in zip(self.position, other.position)
        ))


@dataclass
class Pharmacophore:
    """
    Represents a complete pharmacophore model.

    Attributes:
        features: List of pharmacophore features
        name: Identifier for this pharmacophore
        source_structure: Path to source structure file
    """
    features: List[PharmacophoreFeature] = field(default_factory=list)
    name: str = ""
    source_structure: str = ""

    def add_feature(self, feature: PharmacophoreFeature):
        """Add a feature to the pharmacophore."""
        self.features.append(feature)

    def get_features_by_type(self, feature_type: FeatureType) -> List[PharmacophoreFeature]:
        """Get all features of a specific type."""
        return [f for f in self.features if f.feature_type == feature_type]

    def get_centroid(self) -> Tuple[float, float, float]:
        """Calculate the geometric centroid of all features."""
        if not self.features:
            return (0.0, 0.0, 0.0)

        positions = np.array([f.position for f in self.features])
        centroid = positions.mean(axis=0)
        return tuple(centroid)

    def get_radius(self) -> float:
        """Calculate the radius (max distance from centroid)."""
        if not self.features:
            return 0.0

        centroid = np.array(self.get_centroid())
        positions = np.array([f.position for f in self.features])
        distances = np.linalg.norm(positions - centroid, axis=1)
        return float(distances.max())

    def to_vector(self) -> np.ndarray:
        """Convert pharmacophore to feature vector."""
        return pharmacophore_to_vector(self)


class PharmacophoreGenerator:
    """
    Generates pharmacophore models from molecular structures and interfaces.

    This class can create pharmacophores from:
    - Molecular glue binding interfaces
    - Small molecule ligands
    - Protein surface patches

    Example:
        >>> generator = PharmacophoreGenerator()
        >>> pharmacophore = generator.from_interface_contacts(contacts)
        >>> vector = pharmacophore.to_vector()
    """

    def __init__(self, feature_radius: float = 1.5):
        """
        Initialize the pharmacophore generator.

        Args:
            feature_radius: Default tolerance radius for features (Å)
        """
        self.feature_radius = feature_radius
        self._feat_factory = None

        if Chem is not None and FDEF_FILE is not None:
            try:
                self._feat_factory = ChemicalFeatures.BuildFeatureFactory(FDEF_FILE)
            except Exception:
                pass

    def from_interface_contacts(
        self,
        contacts: List[Dict],
        positions: Optional[Dict[int, Tuple[float, float, float]]] = None
    ) -> Pharmacophore:
        """
        Generate pharmacophore from interface contact residues.

        Args:
            contacts: List of contact residue dictionaries
            positions: Optional mapping of residue IDs to 3D positions

        Returns:
            Pharmacophore object representing the interface
        """
        pharmacophore = Pharmacophore(name="interface_pharmacophore")

        for contact in contacts:
            resname = contact.get('resname', '')
            resid = contact.get('resid', 0)
            distance = contact.get('min_distance', 5.0)

            # Get position (use provided or estimate from distance)
            if positions and resid in positions:
                position = positions[resid]
            else:
                # Placeholder - in real usage, extract from structure
                position = (0.0, 0.0, distance)

            # Weight by distance (closer = more important)
            weight = 1.0 / (1.0 + distance / 5.0)

            # Create features based on residue type
            feature_types = RESIDUE_FEATURES.get(resname, [])

            for feat_type in feature_types:
                feature = PharmacophoreFeature(
                    feature_type=feat_type,
                    position=position,
                    radius=self.feature_radius,
                    source=f"{resname}_{resid}",
                    weight=weight
                )
                pharmacophore.add_feature(feature)

        return pharmacophore

    def from_molecule(self, mol: Union[str, 'Chem.Mol']) -> Pharmacophore:
        """
        Generate pharmacophore from a small molecule.

        Args:
            mol: RDKit molecule object or SMILES/file path

        Returns:
            Pharmacophore object for the molecule
        """
        if Chem is None:
            raise ImportError("RDKit is required for molecule pharmacophores")

        # Handle different input types
        if isinstance(mol, str):
            if mol.endswith('.sdf') or mol.endswith('.mol'):
                mol = Chem.MolFromMolFile(mol)
            elif mol.endswith('.mol2'):
                mol = Chem.MolFromMol2File(mol)
            else:
                mol = Chem.MolFromSmiles(mol)
                if mol is not None:
                    mol = Chem.AddHs(mol)
                    AllChem.EmbedMolecule(mol, randomSeed=42)

        if mol is None:
            raise ValueError("Could not parse molecule")

        pharmacophore = Pharmacophore(name="molecule_pharmacophore")

        # Use RDKit feature factory if available
        if self._feat_factory is not None:
            try:
                feats = self._feat_factory.GetFeaturesForMol(mol)
                conf = mol.GetConformer()

                for feat in feats:
                    family = feat.GetFamily()
                    pos = feat.GetPos()

                    # Map RDKit feature families to our types
                    type_map = {
                        'Donor': FeatureType.HBOND_DONOR,
                        'Acceptor': FeatureType.HBOND_ACCEPTOR,
                        'Aromatic': FeatureType.AROMATIC,
                        'Hydrophobe': FeatureType.HYDROPHOBIC,
                        'PosIonizable': FeatureType.POSITIVE,
                        'NegIonizable': FeatureType.NEGATIVE,
                    }

                    if family in type_map:
                        pharma_feat = PharmacophoreFeature(
                            feature_type=type_map[family],
                            position=(pos.x, pos.y, pos.z),
                            radius=self.feature_radius,
                            source=family,
                            weight=1.0
                        )
                        pharmacophore.add_feature(pharma_feat)

                return pharmacophore

            except Exception as e:
                warnings.warn(f"RDKit pharmacophore failed: {e}, using fallback")

        # Fallback: simple atom-based features
        conf = mol.GetConformer()

        for atom in mol.GetAtoms():
            pos = conf.GetAtomPosition(atom.GetIdx())
            symbol = atom.GetSymbol()

            features = []
            if symbol == 'N' and atom.GetTotalNumHs() > 0:
                features.append(FeatureType.HBOND_DONOR)
            if symbol in ['N', 'O']:
                features.append(FeatureType.HBOND_ACCEPTOR)
            if symbol == 'C' and atom.GetIsAromatic():
                features.append(FeatureType.AROMATIC)
            if symbol == 'C' and not atom.GetIsAromatic():
                features.append(FeatureType.HYDROPHOBIC)

            for feat_type in features:
                pharma_feat = PharmacophoreFeature(
                    feature_type=feat_type,
                    position=(pos.x, pos.y, pos.z),
                    radius=self.feature_radius,
                    source=f"{symbol}_{atom.GetIdx()}",
                    weight=1.0
                )
                pharmacophore.add_feature(pharma_feat)

        return pharmacophore

    def from_protein_surface(
        self,
        surface_residues: List[Dict],
        structure_path: Optional[str] = None
    ) -> Pharmacophore:
        """
        Generate pharmacophore from protein surface residues.

        Args:
            surface_residues: List of surface residue dictionaries
            structure_path: Optional path to structure for coordinates

        Returns:
            Pharmacophore representing the protein surface
        """
        positions = {}

        if structure_path is not None:
            # Extract residue positions from structure
            try:
                import MDAnalysis as mda
                u = mda.Universe(structure_path)
                for res in u.select_atoms("protein").residues:
                    # Use CA position as residue position
                    ca = res.atoms.select_atoms("name CA")
                    if len(ca) > 0:
                        positions[res.resid] = tuple(ca.positions[0])
            except Exception:
                pass

        return self.from_interface_contacts(surface_residues, positions)


def pharmacophore_to_vector(pharmacophore: Pharmacophore) -> np.ndarray:
    """
    Convert a pharmacophore to a fixed-length feature vector.

    The vector encodes:
    - Feature counts by type (7 dimensions)
    - Weighted feature counts (7 dimensions)
    - Spatial distribution statistics (9 dimensions)

    Args:
        pharmacophore: Pharmacophore object to convert

    Returns:
        Numpy array of shape (23,) representing the pharmacophore
    """
    # Initialize feature counts
    feature_counts = {ft: 0 for ft in FeatureType}
    weighted_counts = {ft: 0.0 for ft in FeatureType}

    positions = []
    weights = []

    for feature in pharmacophore.features:
        feature_counts[feature.feature_type] += 1
        weighted_counts[feature.feature_type] += feature.weight
        positions.append(feature.position)
        weights.append(feature.weight)

    # Build vector
    vector = []

    # Feature counts (7 dimensions)
    for ft in FeatureType:
        vector.append(feature_counts[ft])

    # Weighted counts (7 dimensions)
    for ft in FeatureType:
        vector.append(weighted_counts[ft])

    # Spatial statistics (9 dimensions)
    if positions:
        positions = np.array(positions)
        centroid = positions.mean(axis=0)

        # Centroid position
        vector.extend(centroid.tolist())

        # Spread (std in each dimension)
        spread = positions.std(axis=0) if len(positions) > 1 else np.zeros(3)
        vector.extend(spread.tolist())

        # Max extent
        if len(positions) > 1:
            pairwise = euclidean_distances(positions) if euclidean_distances else np.zeros((1, 1))
            max_dist = pairwise.max()
            mean_dist = pairwise.mean()
            volume = np.prod(positions.max(axis=0) - positions.min(axis=0))
        else:
            max_dist = 0.0
            mean_dist = 0.0
            volume = 0.0

        vector.extend([max_dist, mean_dist, volume])
    else:
        # No features - zero padding
        vector.extend([0.0] * 9)

    return np.array(vector, dtype=np.float32)


def calculate_pharmacophore_similarity(
    pharm1: Union[Pharmacophore, np.ndarray],
    pharm2: Union[Pharmacophore, np.ndarray],
    method: str = "cosine"
) -> float:
    """
    Calculate similarity between two pharmacophores.

    Args:
        pharm1: First pharmacophore or its vector representation
        pharm2: Second pharmacophore or its vector representation
        method: Similarity method ("cosine", "euclidean", "tanimoto")

    Returns:
        Similarity score (0-1 for cosine/tanimoto, inverse distance for euclidean)
    """
    # Convert to vectors if needed
    if isinstance(pharm1, Pharmacophore):
        v1 = pharmacophore_to_vector(pharm1)
    else:
        v1 = np.array(pharm1)

    if isinstance(pharm2, Pharmacophore):
        v2 = pharmacophore_to_vector(pharm2)
    else:
        v2 = np.array(pharm2)

    # Calculate similarity
    if method == "cosine":
        if cosine_similarity is not None:
            return float(cosine_similarity(
                v1.reshape(1, -1),
                v2.reshape(1, -1)
            )[0, 0])
        else:
            # Manual cosine similarity
            dot = np.dot(v1, v2)
            norm1 = np.linalg.norm(v1)
            norm2 = np.linalg.norm(v2)
            if norm1 == 0 or norm2 == 0:
                return 0.0
            return float(dot / (norm1 * norm2))

    elif method == "euclidean":
        dist = np.linalg.norm(v1 - v2)
        return float(1.0 / (1.0 + dist))

    elif method == "tanimoto":
        # Tanimoto for count vectors
        intersection = np.minimum(v1, v2).sum()
        union = np.maximum(v1, v2).sum()
        if union == 0:
            return 0.0
        return float(intersection / union)

    else:
        raise ValueError(f"Unknown similarity method: {method}")


def align_pharmacophores(
    reference: Pharmacophore,
    mobile: Pharmacophore,
    feature_types: Optional[List[FeatureType]] = None
) -> Tuple[np.ndarray, float]:
    """
    Align two pharmacophores and return transformation matrix.

    Args:
        reference: Reference pharmacophore (fixed)
        mobile: Mobile pharmacophore (to be aligned)
        feature_types: Optional list of feature types to use for alignment

    Returns:
        Tuple of (4x4 transformation matrix, RMSD after alignment)
    """
    # Filter features by type if specified
    if feature_types:
        ref_feats = [f for f in reference.features if f.feature_type in feature_types]
        mob_feats = [f for f in mobile.features if f.feature_type in feature_types]
    else:
        ref_feats = reference.features
        mob_feats = mobile.features

    if len(ref_feats) < 3 or len(mob_feats) < 3:
        warnings.warn("Not enough features for alignment (need >= 3)")
        return np.eye(4), float('inf')

    # Get positions
    ref_pos = np.array([f.position for f in ref_feats])
    mob_pos = np.array([f.position for f in mob_feats])

    # Use first 3 features for simple alignment
    # More sophisticated approaches would use all features with optimization
    n = min(len(ref_pos), len(mob_pos), 10)

    # Calculate centroids
    ref_centroid = ref_pos[:n].mean(axis=0)
    mob_centroid = mob_pos[:n].mean(axis=0)

    # Center coordinates
    ref_centered = ref_pos[:n] - ref_centroid
    mob_centered = mob_pos[:n] - mob_centroid

    # SVD for rotation
    H = mob_centered.T @ ref_centered
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T

    # Handle reflection
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T

    # Calculate translation
    t = ref_centroid - R @ mob_centroid

    # Build transformation matrix
    transform = np.eye(4)
    transform[:3, :3] = R
    transform[:3, 3] = t

    # Calculate RMSD
    mob_transformed = (R @ mob_pos[:n].T).T + t
    rmsd = np.sqrt(((ref_pos[:n] - mob_transformed) ** 2).mean())

    return transform, float(rmsd)
