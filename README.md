![image](https://github.com/user-attachments/assets/3ef7d5c2-30ff-42f4-9516-6e8ed80c510f)


## Prepare folders for BTK project
`mkdir -p BTK/{ligprep,proteinpreparation,docking,analysis,docs,fep,DNN}`

## STEP 1: Prepare ligands for docking in PDBQT format;  generate 3D conformation using RDKiT, assign protonation states using scrub.py; 

`conda activate autodock`

`cd ligprep`

`bash ligprep.sh` 


## STEP 2: Prepare protein for docking in PDBQT format; fill missing atoms and residues, protein preparation using PDBFixer, OPENMM, PDB2PQR

`protein_preparation_pdbfixer.ipynb`

`bash protein_preparation_pdb_to_pdbqt.sh`

## STEP 3: Docking using Autodock Vina program

`cd ../docking`

`bash vina_docking.sh`

`bash split_results.sh`

## STEP 4: Analysis; select the best docking pose (loweest docking score kcal/mol), visualize the protein-ligand interaction fingerprints 
autodock_results_view.ipynb will combine the protein and ligand to a single PDB

`cd ../analysis`

`bash autodock_results.sh`

## STEP 5: Check the docked poses, sterocenters using Posebusters
Save the bioactive conformation of ligand from the PDB that will be used as reference to calculate the RMSD.

`cd ../analysis`

`bash posebuster.sh`

## STEP 6: FEP Calculation; Docking provides the best pose of the BTK design; FEP calculation prioritize these deisgns

`cd ../fep`

## STEP 7: Deep Neural Network to predict the binding affinity values of designs

`cd ../DNN`

`btk_pic50_DNN.ipynb`


