# Molecular docking steps for BTK molecular glues

# Prepare folders for BTK project
mkdir -p BTK/{ligprep,proteinpreparation,docking,analysis,docs,fep}

# Step1: Ligprep using MEEKO, scrup.py

conda activate autodock

cd ligprep

bash ligprep.sh 


# Step2: Protein preparation using OPENMM, PDBFixer, PDB2PQR

protein_preparation_pdbfixer.ipynb

protein_preparation_pdb_to_pdbqt.sh

# Step3: Docking using Vina

cd ../docking

bash vina_docking.sh

# Step4: Analysis, autodock_results_view.ipynb will combine the protein and ligand to a single PDB

cd ../analysis

bash autodock_results.sh

# Step5: Posebusters

#load the PDB structure in pymol and save the lignd only pdb file.

cd ../analysis

bash posebuster.sh
