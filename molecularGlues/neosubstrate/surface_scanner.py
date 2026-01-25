"""
Surface Scanner Module for Neo-Substrate Discovery

This module provides tools for scanning protein surfaces to identify
patches that are compatible with molecular glue binding. It uses a
computational proteomics approach to compare surface features against
known binding interfaces.

Key Components:
    - ProteinSurface: Represents a protein's solvent-accessible surface
    - SurfaceScanner: Scans proteins for compatible binding patches
    - Proteome scanning utilities for large-scale screening
"""

import numpy as np
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional, Iterator, Set
from pathlib import Path
import warnings
import json
from concurrent.futures import ThreadPoolExecutor, as_completed

try:
    import MDAnalysis as mda
    from MDAnalysis.analysis import distances
except ImportError:
    mda = None
    warnings.warn("MDAnalysis not installed. Surface analysis limited.")

try:
    from scipy.spatial import ConvexHull, Delaunay
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import pdist
except ImportError:
    ConvexHull = None
    warnings.warn("SciPy not installed. Surface calculations limited.")

from .pharmacophore import (
    Pharmacophore,
    PharmacophoreFeature,
    PharmacophoreGenerator,
    pharmacophore_to_vector,
    calculate_pharmacophore_similarity,
    FeatureType,
    RESIDUE_FEATURES,
)


@dataclass
class SurfaceResidue:
    """Represents a residue exposed on the protein surface."""
    resid: int
    resname: str
    chain: str
    position: Tuple[float, float, float]
    sasa: float = 0.0  # Solvent accessible surface area
    features: List[FeatureType] = field(default_factory=list)

    def __post_init__(self):
        """Assign features based on residue type."""
        if not self.features:
            self.features = list(RESIDUE_FEATURES.get(self.resname, []))


@dataclass
class SurfacePatch:
    """
    Represents a contiguous patch on a protein surface.

    A surface patch is a cluster of nearby residues that could potentially
    serve as a binding site for molecular glues or Neo-substrate interactions.
    """
    patch_id: int
    residues: List[SurfaceResidue]
    centroid: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    area: float = 0.0
    pharmacophore: Optional[Pharmacophore] = None

    def __post_init__(self):
        """Calculate derived properties."""
        if self.residues and self.centroid == (0.0, 0.0, 0.0):
            positions = np.array([r.position for r in self.residues])
            self.centroid = tuple(positions.mean(axis=0))

        if self.residues and self.area == 0.0:
            self.area = sum(r.sasa for r in self.residues)

    def get_feature_counts(self) -> Dict[FeatureType, int]:
        """Count features in this patch."""
        counts = {ft: 0 for ft in FeatureType}
        for res in self.residues:
            for feat in res.features:
                counts[feat] += 1
        return counts

    def to_pharmacophore(self) -> Pharmacophore:
        """Convert patch to a pharmacophore model."""
        if self.pharmacophore is not None:
            return self.pharmacophore

        pharmacophore = Pharmacophore(name=f"patch_{self.patch_id}")

        for res in self.residues:
            for feat_type in res.features:
                feature = PharmacophoreFeature(
                    feature_type=feat_type,
                    position=res.position,
                    radius=1.5,
                    source=f"{res.resname}_{res.resid}",
                    weight=1.0
                )
                pharmacophore.add_feature(feature)

        self.pharmacophore = pharmacophore
        return pharmacophore


@dataclass
class ProteinSurface:
    """
    Represents the solvent-accessible surface of a protein.

    Attributes:
        protein_id: Identifier for this protein
        structure_path: Path to structure file
        patches: List of identified surface patches
        total_sasa: Total solvent accessible surface area
    """
    protein_id: str
    structure_path: str = ""
    patches: List[SurfacePatch] = field(default_factory=list)
    total_sasa: float = 0.0
    surface_residues: List[SurfaceResidue] = field(default_factory=list)

    def get_patch_by_residues(self, resids: Set[int]) -> Optional[SurfacePatch]:
        """Find patch containing specified residues."""
        for patch in self.patches:
            patch_resids = {r.resid for r in patch.residues}
            if resids.issubset(patch_resids):
                return patch
        return None


class SurfaceScanner:
    """
    Scans protein surfaces for Neo-substrate compatible patches.

    This class identifies surface patches on proteins that match the
    pharmacophore requirements of a molecular glue interface.

    Example:
        >>> scanner = SurfaceScanner()
        >>> surface = scanner.analyze_protein("protein.pdb")
        >>> compatible = scanner.find_compatible_patches(
        ...     surface,
        ...     reference_pharmacophore,
        ...     similarity_threshold=0.7
        ... )
    """

    def __init__(
        self,
        sasa_threshold: float = 10.0,
        patch_distance: float = 8.0,
        min_patch_size: int = 3,
        max_patch_size: int = 20,
    ):
        """
        Initialize the surface scanner.

        Args:
            sasa_threshold: Minimum SASA for a residue to be considered surface (Å²)
            patch_distance: Maximum distance for residues in same patch (Å)
            min_patch_size: Minimum number of residues in a patch
            max_patch_size: Maximum number of residues in a patch
        """
        self.sasa_threshold = sasa_threshold
        self.patch_distance = patch_distance
        self.min_patch_size = min_patch_size
        self.max_patch_size = max_patch_size
        self._pharmacophore_gen = PharmacophoreGenerator()

    def _calculate_sasa(self, universe: 'mda.Universe') -> Dict[int, float]:
        """
        Calculate per-residue SASA (approximate).

        Uses atom count and exposure as a simple proxy for SASA.
        """
        sasa_values = {}
        protein = universe.select_atoms("protein")

        for residue in protein.residues:
            # Simple approximation: count atoms * average atom surface
            # More accurate would use MDAnalysis SASA or FreeSASA
            n_atoms = len(residue.atoms)

            # Check how many atoms are on surface (not buried)
            surface_atoms = 0
            for atom in residue.atoms:
                # Count nearby atoms as metric of burial
                nearby = universe.select_atoms(
                    f"protein and not resid {residue.resid} and around 4 index {atom.index}"
                )
                if len(nearby) < 15:  # Less buried
                    surface_atoms += 1

            # Approximate SASA based on surface atom count
            sasa_values[residue.resid] = surface_atoms * 15.0  # ~15 Å² per surface atom

        return sasa_values

    def identify_surface_residues(
        self,
        structure_path: str,
        chain: Optional[str] = None
    ) -> List[SurfaceResidue]:
        """
        Identify residues on the protein surface.

        Args:
            structure_path: Path to PDB file
            chain: Optional chain ID to analyze

        Returns:
            List of SurfaceResidue objects
        """
        if mda is None:
            raise ImportError("MDAnalysis is required for surface analysis")

        universe = mda.Universe(structure_path)

        if chain:
            selection = f"protein and (segid {chain} or chainID {chain})"
        else:
            selection = "protein"

        protein = universe.select_atoms(selection)

        # Calculate SASA
        sasa_values = self._calculate_sasa(universe)

        surface_residues = []

        for residue in protein.residues:
            sasa = sasa_values.get(residue.resid, 0.0)

            if sasa >= self.sasa_threshold:
                # Get CA position as residue center
                ca = residue.atoms.select_atoms("name CA")
                if len(ca) > 0:
                    position = tuple(ca.positions[0])
                else:
                    position = tuple(residue.atoms.positions.mean(axis=0))

                surface_res = SurfaceResidue(
                    resid=residue.resid,
                    resname=residue.resname,
                    chain=residue.segid if residue.segid else "A",
                    position=position,
                    sasa=sasa
                )
                surface_residues.append(surface_res)

        return surface_residues

    def cluster_into_patches(
        self,
        surface_residues: List[SurfaceResidue]
    ) -> List[SurfacePatch]:
        """
        Cluster surface residues into contiguous patches.

        Args:
            surface_residues: List of surface residues to cluster

        Returns:
            List of SurfacePatch objects
        """
        if len(surface_residues) < self.min_patch_size:
            return []

        # Get positions
        positions = np.array([r.position for r in surface_residues])

        # Hierarchical clustering
        if ConvexHull is not None and len(positions) > 1:
            try:
                dist_matrix = pdist(positions)
                Z = linkage(dist_matrix, method='average')
                clusters = fcluster(Z, self.patch_distance, criterion='distance')
            except Exception:
                # Fallback to simple distance-based clustering
                clusters = self._simple_cluster(positions)
        else:
            clusters = self._simple_cluster(positions)

        # Group residues by cluster
        patches = []
        unique_clusters = set(clusters)

        for cluster_id in unique_clusters:
            cluster_residues = [
                r for r, c in zip(surface_residues, clusters)
                if c == cluster_id
            ]

            if self.min_patch_size <= len(cluster_residues) <= self.max_patch_size:
                patch = SurfacePatch(
                    patch_id=int(cluster_id),
                    residues=cluster_residues
                )
                patches.append(patch)

        return patches

    def _simple_cluster(self, positions: np.ndarray) -> np.ndarray:
        """Simple distance-based clustering fallback."""
        n = len(positions)
        clusters = np.zeros(n, dtype=int)
        current_cluster = 1

        for i in range(n):
            if clusters[i] == 0:
                clusters[i] = current_cluster

                # Find all points within distance
                for j in range(i + 1, n):
                    if clusters[j] == 0:
                        dist = np.linalg.norm(positions[i] - positions[j])
                        if dist < self.patch_distance:
                            clusters[j] = current_cluster

                current_cluster += 1

        return clusters

    def analyze_protein(
        self,
        structure_path: str,
        protein_id: Optional[str] = None,
        chain: Optional[str] = None
    ) -> ProteinSurface:
        """
        Analyze a protein structure and identify surface patches.

        Args:
            structure_path: Path to PDB file
            protein_id: Identifier for this protein
            chain: Optional chain to analyze

        Returns:
            ProteinSurface object with identified patches
        """
        if protein_id is None:
            protein_id = Path(structure_path).stem

        # Identify surface residues
        surface_residues = self.identify_surface_residues(structure_path, chain)

        # Cluster into patches
        patches = self.cluster_into_patches(surface_residues)

        # Generate pharmacophores for each patch
        for patch in patches:
            patch.to_pharmacophore()

        # Calculate total SASA
        total_sasa = sum(r.sasa for r in surface_residues)

        return ProteinSurface(
            protein_id=protein_id,
            structure_path=structure_path,
            patches=patches,
            total_sasa=total_sasa,
            surface_residues=surface_residues
        )

    def find_compatible_patches(
        self,
        protein_surface: ProteinSurface,
        reference_pharmacophore: Pharmacophore,
        similarity_threshold: float = 0.5,
        method: str = "cosine"
    ) -> List[Tuple[SurfacePatch, float]]:
        """
        Find surface patches compatible with a reference pharmacophore.

        Args:
            protein_surface: ProteinSurface to scan
            reference_pharmacophore: Pharmacophore to match against
            similarity_threshold: Minimum similarity score
            method: Similarity method ("cosine", "tanimoto", "euclidean")

        Returns:
            List of (patch, similarity_score) tuples, sorted by score
        """
        ref_vector = pharmacophore_to_vector(reference_pharmacophore)
        compatible = []

        for patch in protein_surface.patches:
            patch_pharmacophore = patch.to_pharmacophore()
            patch_vector = pharmacophore_to_vector(patch_pharmacophore)

            similarity = calculate_pharmacophore_similarity(
                ref_vector,
                patch_vector,
                method=method
            )

            if similarity >= similarity_threshold:
                compatible.append((patch, similarity))

        # Sort by similarity (highest first)
        compatible.sort(key=lambda x: x[1], reverse=True)

        return compatible


@dataclass
class NeoSubstrateHit:
    """Represents a potential Neo-substrate candidate."""
    protein_id: str
    patch: SurfacePatch
    similarity_score: float
    structure_path: str = ""
    rank: int = 0


def scan_proteome_for_neo_substrates(
    glue_pharmacophore: Pharmacophore,
    proteome_structures: List[str],
    similarity_threshold: float = 0.5,
    n_workers: int = 4,
    progress_callback: Optional[callable] = None
) -> List[NeoSubstrateHit]:
    """
    Scan a proteome for compatible Neo-substrate candidates.

    This is the main entry point for proteome-wide Neo-substrate screening.
    It parallelizes the analysis across multiple protein structures.

    Args:
        glue_pharmacophore: Reference pharmacophore from molecular glue interface
        proteome_structures: List of paths to protein structure files
        similarity_threshold: Minimum similarity for hits
        n_workers: Number of parallel workers
        progress_callback: Optional callback(current, total) for progress

    Returns:
        List of NeoSubstrateHit objects, sorted by similarity score

    Example:
        >>> from pathlib import Path
        >>> structures = list(Path("proteome/").glob("*.pdb"))
        >>> hits = scan_proteome_for_neo_substrates(
        ...     glue_pharmacophore=reference_pharm,
        ...     proteome_structures=[str(s) for s in structures],
        ...     similarity_threshold=0.6
        ... )
        >>> print(f"Found {len(hits)} potential Neo-substrates")
    """
    scanner = SurfaceScanner()
    hits = []
    total = len(proteome_structures)

    def analyze_one(structure_path: str) -> List[NeoSubstrateHit]:
        """Analyze a single protein structure."""
        try:
            surface = scanner.analyze_protein(structure_path)
            compatible = scanner.find_compatible_patches(
                surface,
                glue_pharmacophore,
                similarity_threshold=similarity_threshold
            )

            return [
                NeoSubstrateHit(
                    protein_id=surface.protein_id,
                    patch=patch,
                    similarity_score=score,
                    structure_path=structure_path
                )
                for patch, score in compatible
            ]
        except Exception as e:
            warnings.warn(f"Failed to analyze {structure_path}: {e}")
            return []

    # Parallel processing
    if n_workers > 1:
        with ThreadPoolExecutor(max_workers=n_workers) as executor:
            futures = {
                executor.submit(analyze_one, path): path
                for path in proteome_structures
            }

            for i, future in enumerate(as_completed(futures)):
                result = future.result()
                hits.extend(result)

                if progress_callback:
                    progress_callback(i + 1, total)
    else:
        for i, path in enumerate(proteome_structures):
            result = analyze_one(path)
            hits.extend(result)

            if progress_callback:
                progress_callback(i + 1, total)

    # Sort by similarity and assign ranks
    hits.sort(key=lambda x: x.similarity_score, reverse=True)
    for i, hit in enumerate(hits):
        hit.rank = i + 1

    return hits


def load_proteome_database(database_path: str) -> Dict:
    """
    Load a pre-computed proteome surface database.

    Args:
        database_path: Path to JSON database file

    Returns:
        Dictionary with protein surfaces and metadata
    """
    with open(database_path, 'r') as f:
        data = json.load(f)

    return data


def save_proteome_database(
    surfaces: List[ProteinSurface],
    database_path: str
):
    """
    Save computed protein surfaces to a database file.

    Args:
        surfaces: List of ProteinSurface objects
        database_path: Output path for JSON database
    """
    data = {
        'version': '1.0',
        'n_proteins': len(surfaces),
        'proteins': []
    }

    for surface in surfaces:
        protein_data = {
            'protein_id': surface.protein_id,
            'structure_path': surface.structure_path,
            'total_sasa': surface.total_sasa,
            'n_patches': len(surface.patches),
            'patches': []
        }

        for patch in surface.patches:
            patch_data = {
                'patch_id': patch.patch_id,
                'centroid': list(patch.centroid),
                'area': patch.area,
                'n_residues': len(patch.residues),
                'residues': [
                    {
                        'resid': r.resid,
                        'resname': r.resname,
                        'chain': r.chain,
                        'position': list(r.position),
                        'sasa': r.sasa
                    }
                    for r in patch.residues
                ],
                'pharmacophore_vector': pharmacophore_to_vector(
                    patch.to_pharmacophore()
                ).tolist()
            }
            protein_data['patches'].append(patch_data)

        data['proteins'].append(protein_data)

    with open(database_path, 'w') as f:
        json.dump(data, f, indent=2)
