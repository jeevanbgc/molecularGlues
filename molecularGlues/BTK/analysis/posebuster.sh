#!/bin/bash

for sdf_file in *.sdf; do
    output_file="${sdf_file%.sdf}_bust_results.csv"
    bust "$sdf_file" -l ../proteinpreparation/PDB/ligand_9AJ.pdb -p ../proteinpreparation/protein_prepared/prepared_5vfi.pdb --outfmt csv --output "$output_file"
done
