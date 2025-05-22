#!/bin/bash


output_folder="docking_results"
mkdir -p "$output_folder"

# We used getBox plugin in pymol to get the center and size of the box
for ligand in ../ligprep/mols_pdbqt/*.pdbqt; do
    output_file="$output_folder/$(basename "$ligand" .pdbqt)_out.pdbqt"
    vina --receptor ../proteinpreparation/prepared_5vfi.pdbqt --ligand "$ligand" --center_x 22.4 --center_y 7.9 --center_z 1.3 --size_x 20.1 --size_y 20.6 --size_z 31.2 --seed -26590 --out "$output_file"
done
