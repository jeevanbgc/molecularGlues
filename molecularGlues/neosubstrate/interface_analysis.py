"""
Interface Analysis Module for Molecular Glue Ternary Complexes

This module provides tools for analyzing the interface between molecular glues,
E3 ligases, and potential Neo-substrates. It identifies key contact residues,
calculates interface areas, and characterizes the binding mode.

Key Functions:
    - identify_glue_interface: Find residues contacting molecular glue
    - calculate_interface_area: Compute buried surface area
    - get_contact_residues: Extract detailed contact information
    - InterfaceAnalyzer: Main class for comprehensive interface analysis
"""

import numpy as np
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional, Set
from pathlib import Path
import warnings

try:
    import MDAnalysis as mda
    from MDAnalysis.analysis import distances
    from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
except ImportError:
    mda = None
    warnings.warn("MDAnalysis not installed. Some features will be unavailable.")

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors
except ImportError:
    Chem = None
    warnings.warn("RDKit not installed. Some features will be unavailable.")

try:
    import prolif as plf
except ImportError:
    plf = None
    warnings.warn("ProLIF not installed. Fingerprinting features will be unavailable.")


# Residue classification for interface analysis
RESIDUE_PROPERTIES = {
    'hydrophobic': {'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TRP', 'PRO'},
    'aromatic': {'PHE', 'TYR', 'TRP', 'HIS'},
    'polar': {'SER', 'THR', 'ASN', 'GLN', 'TYR', 'CYS'},
    'positive': {'LYS', 'ARG', 'HIS'},
    'negative': {'ASP', 'GLU'},
    'hbond_donor': {'SER', 'THR', 'TYR', 'ASN', 'GLN', 'HIS', 'LYS', 'ARG', 'TRP'},
    'hbond_acceptor': {'SER', 'THR', 'TYR', 'ASN', 'GLN', 'HIS', 'ASP', 'GLU'},
    'small': {'GLY', 'ALA', 'SER'},
    'bulky': {'PHE', 'TYR', 'TRP', 'ARG', 'LYS'},
}


@dataclass
class ContactResidue:
    """Represents a residue in contact with the molecular glue."""
    resid: int
    resname: str
    chain: str
    min_distance: float
    contact_atoms: List[str] = field(default_factory=list)
    properties: Set[str] = field(default_factory=set)

    def __post_init__(self):
        """Assign residue properties based on residue name."""
        for prop, residues in RESIDUE_PROPERTIES.items():
            if self.resname in residues:
                self.properties.add(prop)


@dataclass
class InterfaceRegion:
    """Represents a binding interface region."""
    e3_contacts: List[ContactResidue]
    substrate_contacts: List[ContactResidue]
    glue_atoms: List[str]
    buried_surface_area: float = 0.0
    hydrogen_bonds: List[Dict] = field(default_factory=list)
    hydrophobic_contacts: int = 0
    electrostatic_interactions: int = 0


class InterfaceAnalyzer:
    """
    Comprehensive analyzer for molecular glue ternary complex interfaces.

    This class provides methods to:
    - Identify contact residues between molecular glue and proteins
    - Calculate interface properties (area, polarity, shape complementarity)
    - Extract pharmacophore-relevant features
    - Compare interfaces across different complexes

    Attributes:
        structure_path: Path to the ternary complex PDB file
        glue_selection: MDAnalysis selection string for molecular glue
        e3_selection: MDAnalysis selection string for E3 ligase
        substrate_selection: MDAnalysis selection string for substrate
        distance_cutoff: Distance threshold for contact definition (Angstroms)

    Example:
        >>> analyzer = InterfaceAnalyzer(
        ...     structure_path="ternary_complex.pdb",
        ...     glue_selection="resname LIG",
        ...     e3_selection="segid A",
        ...     substrate_selection="segid B"
        ... )
        >>> interface = analyzer.analyze()
        >>> print(f"Interface area: {interface.buried_surface_area:.1f} Å²")
    """

    def __init__(
        self,
        structure_path: str,
        glue_selection: str = "resname LIG",
        e3_selection: str = "segid A or chain A",
        substrate_selection: str = "segid B or chain B",
        distance_cutoff: float = 5.0,
    ):
        """
        Initialize the interface analyzer.

        Args:
            structure_path: Path to PDB/PDBx file of ternary complex
            glue_selection: MDAnalysis selection for molecular glue
            e3_selection: MDAnalysis selection for E3 ligase chain(s)
            substrate_selection: MDAnalysis selection for substrate chain(s)
            distance_cutoff: Maximum distance for contact definition (Å)
        """
        self.structure_path = Path(structure_path)
        self.glue_selection = glue_selection
        self.e3_selection = e3_selection
        self.substrate_selection = substrate_selection
        self.distance_cutoff = distance_cutoff

        self._universe = None
        self._glue = None
        self._e3 = None
        self._substrate = None

    def _load_structure(self):
        """Load structure into MDAnalysis Universe."""
        if mda is None:
            raise ImportError("MDAnalysis is required for interface analysis")

        if self._universe is None:
            self._universe = mda.Universe(str(self.structure_path))
            self._glue = self._universe.select_atoms(self.glue_selection)
            self._e3 = self._universe.select_atoms(self.e3_selection)
            self._substrate = self._universe.select_atoms(self.substrate_selection)

            if len(self._glue) == 0:
                raise ValueError(f"No atoms found for glue selection: {self.glue_selection}")

    def get_contact_residues(
        self,
        protein_selection: str = "protein",
        distance_cutoff: Optional[float] = None
    ) -> List[ContactResidue]:
        """
        Identify residues in contact with the molecular glue.

        Args:
            protein_selection: MDAnalysis selection for protein atoms
            distance_cutoff: Override default distance cutoff

        Returns:
            List of ContactResidue objects for all contacting residues
        """
        self._load_structure()

        cutoff = distance_cutoff or self.distance_cutoff
        protein = self._universe.select_atoms(protein_selection)

        contact_residues = []
        seen_residues = set()

        for residue in protein.residues:
            # Calculate minimum distance between glue and residue
            dist_matrix = distances.distance_array(
                self._glue.positions,
                residue.atoms.positions
            )
            min_dist = dist_matrix.min()

            if min_dist < cutoff:
                reskey = (residue.resid, residue.segid)
                if reskey not in seen_residues:
                    seen_residues.add(reskey)

                    # Find contacting atoms
                    contact_atoms = []
                    for i, atom in enumerate(residue.atoms):
                        if dist_matrix[:, i].min() < cutoff:
                            contact_atoms.append(atom.name)

                    contact = ContactResidue(
                        resid=residue.resid,
                        resname=residue.resname,
                        chain=residue.segid if residue.segid else "A",
                        min_distance=float(min_dist),
                        contact_atoms=contact_atoms
                    )
                    contact_residues.append(contact)

        # Sort by distance
        contact_residues.sort(key=lambda x: x.min_distance)

        return contact_residues

    def get_e3_contacts(self) -> List[ContactResidue]:
        """Get E3 ligase residues contacting the molecular glue."""
        self._load_structure()
        return self.get_contact_residues(self.e3_selection)

    def get_substrate_contacts(self) -> List[ContactResidue]:
        """Get substrate residues contacting the molecular glue."""
        self._load_structure()
        return self.get_contact_residues(self.substrate_selection)

    def calculate_interface_area(self) -> Tuple[float, float, float]:
        """
        Calculate buried surface area at the interface.

        Uses the SASA (Solvent Accessible Surface Area) difference method:
        BSA = SASA(separate) - SASA(complex)

        Returns:
            Tuple of (total_bsa, e3_bsa, substrate_bsa) in Å²
        """
        self._load_structure()

        # Approximate BSA using contact atom count * average atom area
        # More accurate calculation would require FreeSASA or similar
        avg_atom_area = 15.0  # Approximate Å² per atom

        e3_contacts = self.get_e3_contacts()
        substrate_contacts = self.get_substrate_contacts()

        e3_contact_atoms = sum(len(c.contact_atoms) for c in e3_contacts)
        substrate_contact_atoms = sum(len(c.contact_atoms) for c in substrate_contacts)

        e3_bsa = e3_contact_atoms * avg_atom_area
        substrate_bsa = substrate_contact_atoms * avg_atom_area
        total_bsa = e3_bsa + substrate_bsa

        return total_bsa, e3_bsa, substrate_bsa

    def analyze_hydrogen_bonds(self) -> List[Dict]:
        """
        Identify hydrogen bonds at the interface.

        Returns:
            List of hydrogen bond dictionaries with donor, acceptor, distance, angle
        """
        self._load_structure()

        hbonds = []

        # Simple geometric H-bond detection
        # More sophisticated analysis would use HydrogenBondAnalysis
        donors = self._universe.select_atoms(
            f"({self.glue_selection}) and (name N* O*)"
        )
        acceptors = self._universe.select_atoms(
            f"protein and (name O* N*) and around {self.distance_cutoff} ({self.glue_selection})"
        )

        for donor in donors:
            for acceptor in acceptors:
                dist = np.linalg.norm(donor.position - acceptor.position)
                if 2.5 < dist < 3.5:  # Typical H-bond distance
                    hbonds.append({
                        'donor': f"{donor.resname}_{donor.resid}:{donor.name}",
                        'acceptor': f"{acceptor.resname}_{acceptor.resid}:{acceptor.name}",
                        'distance': float(dist),
                        'type': 'potential_hbond'
                    })

        return hbonds

    def analyze(self) -> InterfaceRegion:
        """
        Perform comprehensive interface analysis.

        Returns:
            InterfaceRegion object containing all analysis results
        """
        self._load_structure()

        e3_contacts = self.get_e3_contacts()
        substrate_contacts = self.get_substrate_contacts()
        total_bsa, _, _ = self.calculate_interface_area()
        hbonds = self.analyze_hydrogen_bonds()

        # Count interaction types
        hydrophobic = sum(
            1 for c in e3_contacts + substrate_contacts
            if 'hydrophobic' in c.properties
        )
        electrostatic = sum(
            1 for c in e3_contacts + substrate_contacts
            if 'positive' in c.properties or 'negative' in c.properties
        )

        glue_atoms = [atom.name for atom in self._glue]

        return InterfaceRegion(
            e3_contacts=e3_contacts,
            substrate_contacts=substrate_contacts,
            glue_atoms=glue_atoms,
            buried_surface_area=total_bsa,
            hydrogen_bonds=hbonds,
            hydrophobic_contacts=hydrophobic,
            electrostatic_interactions=electrostatic
        )

    def get_interface_fingerprint(self) -> Optional[np.ndarray]:
        """
        Generate ProLIF interaction fingerprint for the interface.

        Returns:
            Numpy array of interaction fingerprint or None if ProLIF unavailable
        """
        if plf is None:
            warnings.warn("ProLIF not available for fingerprinting")
            return None

        self._load_structure()

        try:
            # Create ProLIF molecule objects
            protein_mol = plf.Molecule.from_mda(self._universe.select_atoms("protein"))
            ligand_mol = plf.Molecule.from_mda(self._glue)

            # Generate fingerprint
            fp = plf.Fingerprint()
            fp.run_from_iterable([ligand_mol], protein_mol)

            return fp.to_numpy()
        except Exception as e:
            warnings.warn(f"Fingerprint generation failed: {e}")
            return None

    def to_dict(self) -> Dict:
        """Export analysis results as dictionary."""
        interface = self.analyze()

        return {
            'structure': str(self.structure_path),
            'e3_contacts': [
                {
                    'resid': c.resid,
                    'resname': c.resname,
                    'chain': c.chain,
                    'distance': c.min_distance,
                    'properties': list(c.properties)
                }
                for c in interface.e3_contacts
            ],
            'substrate_contacts': [
                {
                    'resid': c.resid,
                    'resname': c.resname,
                    'chain': c.chain,
                    'distance': c.min_distance,
                    'properties': list(c.properties)
                }
                for c in interface.substrate_contacts
            ],
            'buried_surface_area': interface.buried_surface_area,
            'hydrogen_bonds': interface.hydrogen_bonds,
            'hydrophobic_contacts': interface.hydrophobic_contacts,
            'electrostatic_interactions': interface.electrostatic_interactions
        }


def identify_glue_interface(
    ternary_pdb: str,
    glue_resname: str = "LIG",
    distance_cutoff: float = 5.0
) -> List[ContactResidue]:
    """
    Convenience function to identify residues at the molecular glue interface.

    Args:
        ternary_pdb: Path to ternary complex PDB file
        glue_resname: Residue name of molecular glue
        distance_cutoff: Maximum distance for contact definition

    Returns:
        List of ContactResidue objects
    """
    analyzer = InterfaceAnalyzer(
        structure_path=ternary_pdb,
        glue_selection=f"resname {glue_resname}",
        distance_cutoff=distance_cutoff
    )
    return analyzer.get_contact_residues()


def calculate_interface_area(ternary_pdb: str, glue_resname: str = "LIG") -> float:
    """
    Calculate total buried surface area at molecular glue interface.

    Args:
        ternary_pdb: Path to ternary complex PDB file
        glue_resname: Residue name of molecular glue

    Returns:
        Buried surface area in Å²
    """
    analyzer = InterfaceAnalyzer(
        structure_path=ternary_pdb,
        glue_selection=f"resname {glue_resname}"
    )
    total_bsa, _, _ = analyzer.calculate_interface_area()
    return total_bsa


def get_contact_residues(
    ternary_pdb: str,
    glue_resname: str = "LIG",
    chain: Optional[str] = None,
    distance_cutoff: float = 5.0
) -> List[Dict]:
    """
    Get detailed contact residue information.

    Args:
        ternary_pdb: Path to ternary complex PDB file
        glue_resname: Residue name of molecular glue
        chain: Optional chain ID to filter contacts
        distance_cutoff: Maximum distance for contact definition

    Returns:
        List of dictionaries with contact residue details
    """
    analyzer = InterfaceAnalyzer(
        structure_path=ternary_pdb,
        glue_selection=f"resname {glue_resname}",
        distance_cutoff=distance_cutoff
    )

    if chain:
        selection = f"protein and (segid {chain} or chain {chain})"
    else:
        selection = "protein"

    contacts = analyzer.get_contact_residues(selection)

    return [
        {
            'resid': c.resid,
            'resname': c.resname,
            'chain': c.chain,
            'min_distance': c.min_distance,
            'contact_atoms': c.contact_atoms,
            'properties': list(c.properties)
        }
        for c in contacts
    ]
