#!/bin/bash

export MGL_ROOT=/mnt/data/software_install/mgltools_x86_64Linux2_1.5.7p1/mgltools_x86_64Linux2_1.5.7

# Use the prepared protein uing PDBfixer and convert to pdbqt
$MGL_ROOT/bin/pythonsh $MGL_ROOT/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py -r protein_prepared/prepared_5vfi.pdb  -o prepared_5vfi.pdbqt