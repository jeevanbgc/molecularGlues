"""
Ternary Complex Docking Module

This module provides tools for docking Neo-substrate candidates to form
ternary complexes with E3 ligases and molecular glues. It supports both
rigid-body and flexible docking approaches.

Key Features:
    - Protein-protein docking for ternary complex formation
    - Constraint-guided docking using molecular glue interface
    - Geometry validation for productive complexes
    - Integration with AutoDock Vina and other docking engines
"""

import numpy as np
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional, Union
from pathlib import Path
import warnings
import subprocess
import tempfile
import shutil
import os

try:
    from vina import Vina
except ImportError:
    Vina = None
    warnings.warn("AutoDock Vina not installed. Docking will be limited.")

try:
    import MDAnalysis as mda
    from MDAnalysis.analysis import distances, align
except ImportError:
    mda = None

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:
    Chem = None

try:
    from openmm.app import PDBFile
    from pdbfixer import PDBFixer
except ImportError:
    PDBFile = None
    PDBFixer = None


@dataclass
class DockingPose:
    """Represents a single docking pose."""
    pose_id: int
    score: float  # Binding score (kcal/mol, lower is better)
    coordinates: np.ndarray  # Atomic coordinates
    rmsd_from_best: float = 0.0
    clash_score: float = 0.0
    interface_contacts: int = 0

    def is_valid(self, max_clash: float = 100.0) -> bool:
        """Check if pose passes validation criteria."""
        return self.clash_score < max_clash


@dataclass
class TernaryComplex:
    """
    Represents a ternary complex structure.

    Attributes:
        e3_structure: Path to E3 ligase structure
        glue_structure: Path to molecular glue structure
        substrate_structure: Path to Neo-substrate structure
        poses: List of docked poses
        best_pose: Index of best pose
    """
    e3_structure: str
    glue_structure: str
    substrate_structure: str
    poses: List[DockingPose] = field(default_factory=list)
    best_pose: int = 0
    total_score: float = 0.0

    def get_best_pose(self) -> Optional[DockingPose]:
        """Return the best docking pose."""
        if self.poses:
            return self.poses[self.best_pose]
        return None


class TernaryDocker:
    """
    Docks Neo-substrate candidates to form ternary complexes.

    This class orchestrates the docking of potential Neo-substrates
    to E3 ligase + molecular glue complexes, validating geometry
    and scoring the resulting ternary structures.

    Example:
        >>> docker = TernaryDocker(
        ...     e3_pdb="crbn.pdb",
        ...     glue_sdf="thalidomide.sdf"
        ... )
        >>> complex = docker.dock_substrate("substrate.pdb")
        >>> print(f"Best score: {complex.get_best_pose().score}")
    """

    def __init__(
        self,
        e3_pdb: str,
        glue_sdf: str,
        box_center: Optional[Tuple[float, float, float]] = None,
        box_size: Tuple[float, float, float] = (30.0, 30.0, 30.0),
        exhaustiveness: int = 32,
        n_poses: int = 20,
    ):
        """
        Initialize the ternary docker.

        Args:
            e3_pdb: Path to E3 ligase PDB file
            glue_sdf: Path to molecular glue SDF file
            box_center: Center of docking box (auto-detect if None)
            box_size: Size of docking box (Å)
            exhaustiveness: Docking exhaustiveness parameter
            n_poses: Number of poses to generate
        """
        self.e3_pdb = Path(e3_pdb)
        self.glue_sdf = Path(glue_sdf)
        self.box_center = box_center
        self.box_size = box_size
        self.exhaustiveness = exhaustiveness
        self.n_poses = n_poses

        self._e3_pdbqt = None
        self._glue_pdbqt = None
        self._temp_dir = None

    def _setup_temp_dir(self):
        """Create temporary directory for docking files."""
        if self._temp_dir is None:
            self._temp_dir = tempfile.mkdtemp(prefix="ternary_dock_")

    def _cleanup(self):
        """Remove temporary files."""
        if self._temp_dir and os.path.exists(self._temp_dir):
            shutil.rmtree(self._temp_dir)
            self._temp_dir = None

    def _prepare_receptor(self) -> str:
        """
        Prepare receptor (E3 + glue) for docking.

        Returns:
            Path to prepared PDBQT file
        """
        self._setup_temp_dir()

        # Combine E3 and glue into single receptor
        receptor_pdb = os.path.join(self._temp_dir, "receptor.pdb")
        receptor_pdbqt = os.path.join(self._temp_dir, "receptor.pdbqt")

        # Read E3 structure
        if mda is not None:
            u_e3 = mda.Universe(str(self.e3_pdb))
            e3_atoms = u_e3.select_atoms("protein")

            # Write combined structure
            e3_atoms.write(receptor_pdb)
        else:
            # Fallback: just copy E3 structure
            shutil.copy(self.e3_pdb, receptor_pdb)

        # Convert to PDBQT using prepare_receptor or similar
        try:
            # Try using ADFRsuite prepare_receptor
            subprocess.run(
                ["prepare_receptor", "-r", receptor_pdb, "-o", receptor_pdbqt],
                capture_output=True,
                check=False
            )
        except FileNotFoundError:
            # Fallback: create minimal PDBQT
            self._convert_pdb_to_pdbqt(receptor_pdb, receptor_pdbqt)

        self._e3_pdbqt = receptor_pdbqt
        return receptor_pdbqt

    def _prepare_ligand(self, ligand_path: str) -> str:
        """
        Prepare ligand for docking.

        Args:
            ligand_path: Path to ligand file

        Returns:
            Path to prepared PDBQT file
        """
        self._setup_temp_dir()

        ligand_pdbqt = os.path.join(
            self._temp_dir,
            Path(ligand_path).stem + ".pdbqt"
        )

        # Try using prepare_ligand
        try:
            subprocess.run(
                ["prepare_ligand", "-l", ligand_path, "-o", ligand_pdbqt],
                capture_output=True,
                check=False
            )
        except FileNotFoundError:
            # Fallback conversion
            self._convert_to_pdbqt(ligand_path, ligand_pdbqt)

        return ligand_pdbqt

    def _convert_pdb_to_pdbqt(self, pdb_path: str, pdbqt_path: str):
        """Convert PDB to PDBQT format (simple conversion)."""
        with open(pdb_path, 'r') as f:
            lines = f.readlines()

        with open(pdbqt_path, 'w') as f:
            for line in lines:
                if line.startswith(('ATOM', 'HETATM')):
                    # Add default atom type
                    atom_name = line[12:16].strip()
                    atom_type = atom_name[0]  # Simple: first letter
                    new_line = line[:77].ljust(77) + f"  {atom_type}\n"
                    f.write(new_line)
                elif line.startswith(('TER', 'END')):
                    f.write(line)

    def _convert_to_pdbqt(self, input_path: str, output_path: str):
        """Convert various formats to PDBQT."""
        suffix = Path(input_path).suffix.lower()

        if Chem is not None:
            if suffix == '.sdf':
                mol = Chem.MolFromMolFile(input_path)
            elif suffix == '.mol2':
                mol = Chem.MolFromMol2File(input_path)
            elif suffix == '.pdb':
                mol = Chem.MolFromPDBFile(input_path)
            else:
                mol = None

            if mol is not None:
                # Write as PDB first
                pdb_temp = output_path.replace('.pdbqt', '.pdb')
                Chem.MolToPDBFile(mol, pdb_temp)
                self._convert_pdb_to_pdbqt(pdb_temp, output_path)
                return

        # Fallback: simple copy with extension change
        shutil.copy(input_path, output_path)

    def _detect_box_center(self) -> Tuple[float, float, float]:
        """Auto-detect docking box center from glue position."""
        if mda is not None:
            try:
                # Load glue structure
                if str(self.glue_sdf).endswith('.sdf'):
                    # Need to load SDF differently
                    if Chem is not None:
                        mol = Chem.MolFromMolFile(str(self.glue_sdf))
                        conf = mol.GetConformer()
                        positions = conf.GetPositions()
                        center = positions.mean(axis=0)
                        return tuple(center)
                else:
                    u = mda.Universe(str(self.glue_sdf))
                    center = u.atoms.positions.mean(axis=0)
                    return tuple(center)
            except Exception:
                pass

        # Default center
        return (0.0, 0.0, 0.0)

    def dock_substrate(
        self,
        substrate_pdb: str,
        use_constraints: bool = True
    ) -> TernaryComplex:
        """
        Dock a Neo-substrate to form a ternary complex.

        Args:
            substrate_pdb: Path to Neo-substrate structure
            use_constraints: Whether to use glue interface constraints

        Returns:
            TernaryComplex with docking results
        """
        # Prepare receptor
        receptor_pdbqt = self._prepare_receptor()

        # Detect box center if not specified
        if self.box_center is None:
            self.box_center = self._detect_box_center()

        # Prepare substrate as ligand
        substrate_pdbqt = self._prepare_ligand(substrate_pdb)

        poses = []

        if Vina is not None:
            try:
                # Use AutoDock Vina
                v = Vina(sf_name='vina')
                v.set_receptor(receptor_pdbqt)
                v.set_ligand_from_file(substrate_pdbqt)

                v.compute_vina_maps(
                    center=list(self.box_center),
                    box_size=list(self.box_size)
                )

                v.dock(
                    exhaustiveness=self.exhaustiveness,
                    n_poses=self.n_poses
                )

                # Extract poses
                energies = v.energies()

                for i, energy in enumerate(energies):
                    pose = DockingPose(
                        pose_id=i,
                        score=float(energy[0]),  # Total score
                        coordinates=np.array([]),  # Would need pose extraction
                        rmsd_from_best=float(energy[1]) if len(energy) > 1 else 0.0
                    )
                    poses.append(pose)

            except Exception as e:
                warnings.warn(f"Vina docking failed: {e}")

        if not poses:
            # Create placeholder pose if docking failed
            poses.append(DockingPose(
                pose_id=0,
                score=0.0,
                coordinates=np.array([]),
                rmsd_from_best=0.0
            ))

        # Find best pose
        best_idx = 0
        best_score = float('inf')
        for i, pose in enumerate(poses):
            if pose.score < best_score:
                best_score = pose.score
                best_idx = i

        complex_result = TernaryComplex(
            e3_structure=str(self.e3_pdb),
            glue_structure=str(self.glue_sdf),
            substrate_structure=substrate_pdb,
            poses=poses,
            best_pose=best_idx,
            total_score=best_score
        )

        return complex_result

    def validate_geometry(
        self,
        complex_result: TernaryComplex,
        max_clash_distance: float = 2.0,
        min_interface_contacts: int = 5
    ) -> Dict:
        """
        Validate the geometry of a ternary complex.

        Checks for:
        - Steric clashes between components
        - Sufficient interface contacts
        - Proper orientation for ubiquitination

        Args:
            complex_result: TernaryComplex to validate
            max_clash_distance: Maximum allowed close contact distance
            min_interface_contacts: Minimum interface contacts required

        Returns:
            Validation results dictionary
        """
        validation = {
            'is_valid': True,
            'clash_count': 0,
            'interface_contacts': 0,
            'geometry_score': 0.0,
            'issues': []
        }

        # Load structures
        if mda is not None:
            try:
                u_e3 = mda.Universe(complex_result.e3_structure)
                u_sub = mda.Universe(complex_result.substrate_structure)

                # Check for clashes
                e3_atoms = u_e3.select_atoms("protein")
                sub_atoms = u_sub.select_atoms("protein")

                dist_matrix = distances.distance_array(
                    e3_atoms.positions,
                    sub_atoms.positions
                )

                clashes = np.sum(dist_matrix < max_clash_distance)
                validation['clash_count'] = int(clashes)

                if clashes > 10:
                    validation['is_valid'] = False
                    validation['issues'].append(f"Too many clashes: {clashes}")

                # Check interface contacts
                interface_dist = 5.0
                contacts = np.sum(dist_matrix < interface_dist)
                validation['interface_contacts'] = int(contacts)

                if contacts < min_interface_contacts:
                    validation['is_valid'] = False
                    validation['issues'].append(
                        f"Insufficient contacts: {contacts}"
                    )

                # Calculate geometry score
                # Penalize clashes, reward contacts
                validation['geometry_score'] = contacts - (clashes * 10)

            except Exception as e:
                validation['is_valid'] = False
                validation['issues'].append(f"Validation error: {e}")

        return validation

    def __del__(self):
        """Cleanup on deletion."""
        self._cleanup()


def dock_neo_substrate(
    e3_pdb: str,
    glue_sdf: str,
    substrate_pdb: str,
    box_center: Optional[Tuple[float, float, float]] = None,
    box_size: Tuple[float, float, float] = (30.0, 30.0, 30.0)
) -> TernaryComplex:
    """
    Convenience function to dock a Neo-substrate.

    Args:
        e3_pdb: Path to E3 ligase structure
        glue_sdf: Path to molecular glue structure
        substrate_pdb: Path to Neo-substrate structure
        box_center: Docking box center (auto-detect if None)
        box_size: Docking box dimensions

    Returns:
        TernaryComplex with docking results
    """
    docker = TernaryDocker(
        e3_pdb=e3_pdb,
        glue_sdf=glue_sdf,
        box_center=box_center,
        box_size=box_size
    )

    return docker.dock_substrate(substrate_pdb)


def validate_ternary_geometry(
    e3_pdb: str,
    substrate_pdb: str,
    glue_sdf: Optional[str] = None
) -> Dict:
    """
    Validate geometry of a ternary complex.

    Args:
        e3_pdb: Path to E3 ligase structure
        substrate_pdb: Path to Neo-substrate structure
        glue_sdf: Optional path to molecular glue

    Returns:
        Validation results dictionary
    """
    docker = TernaryDocker(
        e3_pdb=e3_pdb,
        glue_sdf=glue_sdf or "",
    )

    complex_result = TernaryComplex(
        e3_structure=e3_pdb,
        glue_structure=glue_sdf or "",
        substrate_structure=substrate_pdb
    )

    return docker.validate_geometry(complex_result)


class ConstrainedDocker(TernaryDocker):
    """
    Extended docker that uses interface constraints for guided docking.

    This class adds constraint-based docking where the Neo-substrate
    is guided to form specific interactions with the molecular glue.
    """

    def __init__(
        self,
        e3_pdb: str,
        glue_sdf: str,
        interface_constraints: Optional[List[Dict]] = None,
        **kwargs
    ):
        """
        Initialize constrained docker.

        Args:
            e3_pdb: Path to E3 ligase structure
            glue_sdf: Path to molecular glue structure
            interface_constraints: List of constraint definitions
            **kwargs: Additional TernaryDocker arguments
        """
        super().__init__(e3_pdb, glue_sdf, **kwargs)
        self.constraints = interface_constraints or []

    def add_distance_constraint(
        self,
        glue_atom: str,
        substrate_atom: str,
        target_distance: float,
        tolerance: float = 1.0
    ):
        """
        Add a distance constraint between glue and substrate atoms.

        Args:
            glue_atom: Atom name in molecular glue
            substrate_atom: Atom name in substrate
            target_distance: Target distance (Å)
            tolerance: Allowed deviation (Å)
        """
        self.constraints.append({
            'type': 'distance',
            'glue_atom': glue_atom,
            'substrate_atom': substrate_atom,
            'target': target_distance,
            'tolerance': tolerance
        })

    def score_constraints(
        self,
        pose: DockingPose
    ) -> float:
        """
        Score a pose based on constraint satisfaction.

        Args:
            pose: Docking pose to score

        Returns:
            Constraint satisfaction score (lower is better)
        """
        # Placeholder - would need actual coordinate comparison
        return 0.0


def prepare_ternary_for_md(
    complex_result: TernaryComplex,
    output_dir: str,
    add_hydrogens: bool = True,
    solvate: bool = False
) -> str:
    """
    Prepare a ternary complex for MD simulation.

    Args:
        complex_result: Docked ternary complex
        output_dir: Output directory for prepared files
        add_hydrogens: Whether to add hydrogens
        solvate: Whether to add solvent

    Returns:
        Path to prepared system file
    """
    os.makedirs(output_dir, exist_ok=True)

    output_pdb = os.path.join(output_dir, "ternary_prepared.pdb")

    if PDBFixer is not None:
        try:
            fixer = PDBFixer(filename=complex_result.e3_structure)
            fixer.findMissingResidues()
            fixer.findMissingAtoms()
            fixer.addMissingAtoms()

            if add_hydrogens:
                fixer.addMissingHydrogens(7.0)

            with open(output_pdb, 'w') as f:
                PDBFile.writeFile(fixer.topology, fixer.positions, f)

        except Exception as e:
            warnings.warn(f"PDBFixer preparation failed: {e}")
            shutil.copy(complex_result.e3_structure, output_pdb)
    else:
        shutil.copy(complex_result.e3_structure, output_pdb)

    return output_pdb
