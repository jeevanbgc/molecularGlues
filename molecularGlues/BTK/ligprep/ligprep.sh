# !/bin/bash


# prepare the sdf file
#sdf_file="btk_fenebrutinib_analogs.sdf"
#scrub.py $sdf_file --ph_low 6 --ph_high 8 -o mols_ph6_ph8.sdf

mk_prepare_ligand.py -i mols_ph6_ph8_pymol_selected.sdf --multimol_outdir mols_pdbqt
