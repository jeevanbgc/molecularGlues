# Molecular docking steps for BTK molecular glues

# prepare folders for BTK project
mkdir -p BTK/ligprep BTK/proteinpreparation BTK/docking BTK/analysis BTK/docs
OR mkdir -p BTK/{ligprep,proteinpreparation,docking,analysis,docs,fep}

# Step1: Ligprep using MEEKO, scrup.py
# scrub.py Protonate molecules and add 3D coordinates. Visualize prepared ligands using ligprep_visualize.ipynb, then save selected ligands using pymol
# In order to correct amines,  delete H-atoms and change the formal charge to zero, using 
# pymol command line "alter pk1, formal_charge=0" , where pk1 is the picked atom  using CTRL + mouse right
conda activate autodock
cd ligprep
bash ligprep.sh 


# Step2: Protein preparation using OPENMM, PDBFixer, PDB2PQR
# requires openmm and MGLTools
protein_preparation_pdbfixer.ipynb
protein_preparation_pdb_to_pdbqt.sh

# Step3: Docking using Vina
cd ../docking
bash vina_docking.sh

# Step4: Analysis, autodock_results_view.ipynb will combine the protein and ligand to a single PDB
# using MDanalaysis package, prolif to view protein-ligand interactions
# openbabel to convert sdf to pdb
# Use pymol to view interactions, presets --> small molecule binding, hide all non-polar atoms
cd ../analysis
bash autodock_results.sh

# Posebusters
#load the PDB structure in pymol and save the lignd only pdb file.
cd ../analysis
bash posebuster.sh





Additional NOTES
# Check the protonation state of the protein
# https://github.com/janash/iqb-2024/blob/main/docking_preparation.ipynb
!pdb2pqr --pdb-output=protein_structures/protein_h.pdb --pH=7.4 protein_structures/protein_2zq2.pdb protein_structures/protein_2zq2.pqr --whitespace

# make a directory for pdb files
pdbqt_directory = ## fill in the name of the directory to write PDBQT files to
os.makedirs(pdbqt_directory, exist_ok=True)

u = mda.Universe(f"{protein_directory}/protein_{pdb_id}.pqr")
u.atoms.write(f"{pdbqt_directory}/{pdb_id}.pdbqt")


# scrub.py

    generate 3D coordinates using RDKit's ETKDGv3 and UFF minimization
    enumerate tautomers (aiming at low energy states only)
    enumerate pH corrections
    convert boats to chairs (6-member rings) and enumerate both chair states
    enumerate chiral centers (not implemented right now)
    
# Posebuster

MOL_PRED loaded         . 
Sanitization            .
All atoms connected     . 
Bond lengths            Fail
Bond angles             Fail
Internal steric clash   . 
Aromatic ring flatness  .
Double bond flatness    .
Internal energy         .

# Binding site investigations by janash
https://github.com/janash/iqb-2024/blob/main/extras/binding_site_investigation.ipynb

