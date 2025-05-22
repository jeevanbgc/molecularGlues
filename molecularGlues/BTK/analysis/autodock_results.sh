# !/bin/bash

# Export the results, convert to sdf format and move to analysis folder
# Split_results contain the tope rakned pose
docked_pdbqt="../docking/split_results/*.pdbqt"

for file in $docked_pdbqt; do
    mk_export.py "$file" -s "${file%.pdbqt}_vina_out.sdf"
done


mv ../docking/split_results/*.sdf ../analysis