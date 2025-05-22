mkdir -p {prep,run,images,docs}


# Merge ligands to a singel file using pymol



# Partial charge of ligands
conda create --name  -c conda-forge "espaloma=0.3.2"
conda activate espaloma
pip install espaloma_charge
python esp_neutral.py ligands.sdf  # ligands with ML (Espaloma) partial charges

# Login to aws
aws ec2 start-instances --instance-ids *AWS_KEY_ID*
