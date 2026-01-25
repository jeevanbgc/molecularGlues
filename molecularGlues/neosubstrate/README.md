# Neo-Substrate Discovery Module

A computational pipeline for identifying Neo-substrates that form productive ternary complexes with molecular glues and E3 ligases.

## Overview

This module implements a computational proteomics approach to Neo-substrate discovery:

1. **Interface Analysis** - Characterize molecular glue binding interfaces
2. **Pharmacophore Generation** - Extract key interaction features
3. **Surface Scanning** - Screen proteome for compatible binding patches
4. **Ternary Docking** - Validate candidate binding geometry
5. **Scoring & Ranking** - Prioritize candidates for experimental validation

## Installation

The module requires the `autodock` conda environment:

```bash
conda activate autodock
```

All dependencies (MDAnalysis, RDKit, AutoDock Vina, scikit-learn, etc.) are included in the existing environment.

## Quick Start

```python
from neosubstrate import NeoSubstratePipeline

# Initialize pipeline
pipeline = NeoSubstratePipeline(
    e3_structure="crbn.pdb",
    glue_structure="thalidomide.sdf",
    ternary_structure="ternary_complex.pdb"  # Optional
)

# Run discovery
candidates = pipeline.run(
    proteome_path="proteome/",
    top_k=100,
    dock=True,
    validate=True
)

# Export results
pipeline.export_results("neo_substrates.csv")
```

## Module Structure

```
neosubstrate/
├── __init__.py              # Package exports
├── interface_analysis.py    # Glue-protein interface analysis
├── pharmacophore.py         # Pharmacophore generation and matching
├── surface_scanner.py       # Proteome surface scanning
├── ternary_dock.py          # Ternary complex docking
├── scoring.py               # Candidate scoring and ranking
├── pipeline.py              # Main pipeline orchestration
├── neo_substrate_discovery.ipynb  # Tutorial notebook
└── README.md                # This file
```

## Core Components

### Interface Analysis

```python
from neosubstrate import InterfaceAnalyzer

analyzer = InterfaceAnalyzer(
    structure_path="ternary.pdb",
    glue_selection="resname LIG",
    e3_selection="chain A",
    substrate_selection="chain B"
)

interface = analyzer.analyze()
print(f"Contacts: {len(interface.e3_contacts)}")
print(f"BSA: {interface.buried_surface_area:.1f} Å²")
```

### Pharmacophore Generation

```python
from neosubstrate import PharmacophoreGenerator

gen = PharmacophoreGenerator()

# From interface contacts
pharmacophore = gen.from_interface_contacts(contacts)

# From molecule
pharmacophore = gen.from_molecule("ligand.sdf")

# Convert to vector for comparison
vector = pharmacophore.to_vector()
```

### Surface Scanning

```python
from neosubstrate import SurfaceScanner, scan_proteome_for_neo_substrates

# Single protein
scanner = SurfaceScanner()
surface = scanner.analyze_protein("protein.pdb")
compatible = scanner.find_compatible_patches(surface, reference_pharm)

# Proteome-wide
hits = scan_proteome_for_neo_substrates(
    glue_pharmacophore=reference_pharm,
    proteome_structures=["*.pdb"],
    similarity_threshold=0.5
)
```

### Ternary Docking

```python
from neosubstrate import TernaryDocker

docker = TernaryDocker(
    e3_pdb="crbn.pdb",
    glue_sdf="thalidomide.sdf",
    box_size=(30, 30, 30)
)

result = docker.dock_substrate("substrate.pdb")
print(f"Best score: {result.get_best_pose().score} kcal/mol")
```

### Scoring

```python
from neosubstrate import NeoSubstrateScorer, ScoringWeights

weights = ScoringWeights(
    pharmacophore=0.30,
    docking=0.25,
    interface=0.20,
    druggability=0.15,
    structural=0.10
)

scorer = NeoSubstrateScorer(
    reference_pharmacophore=ref_pharm,
    weights=weights
)

scores = scorer.score_candidates(hits)
ranked = scorer.rank(scores)
```

## Workflow Details

### Step 1: Interface Analysis

The interface analyzer identifies:
- E3 ligase residues contacting the molecular glue
- Substrate residues contacting the molecular glue
- Hydrogen bonds at the interface
- Hydrophobic contacts
- Buried surface area

### Step 2: Pharmacophore Generation

Features extracted include:
- Hydrophobic centers
- Aromatic rings
- H-bond donors/acceptors
- Positive/negative ionizable groups

### Step 3: Surface Scanning

The scanner:
1. Identifies solvent-accessible surface residues
2. Clusters residues into contiguous patches
3. Generates pharmacophores for each patch
4. Compares patches to reference using similarity metrics

### Step 4: Ternary Docking

Docking validates:
- Geometric compatibility
- Binding energy
- Steric clashes
- Interface contacts

### Step 5: Scoring

The scoring function combines:
- Pharmacophore similarity (cosine/Tanimoto)
- Docking score (normalized)
- Interface quality (area, contacts)
- Druggability metrics
- Structural features (lysine accessibility)

## Output Formats

Results can be exported as:
- CSV (tabular data)
- TSV (tab-separated)
- JSON (full details)

## Example Use Cases

### 1. CRBN-based Molecular Glue

```python
pipeline = NeoSubstratePipeline(
    e3_structure="6h0f_crbn.pdb",
    glue_structure="lenalidomide.sdf"
)
candidates = pipeline.run("human_proteome/")
```

### 2. VHL-based Molecular Glue

```python
pipeline = NeoSubstratePipeline(
    e3_structure="3zrc_vhl.pdb",
    glue_structure="vhl_ligand.sdf"
)
candidates = pipeline.run("kinome/")
```

### 3. Quick Screening (No Docking)

```python
from neosubstrate.pipeline import run_quick_scan

hits = run_quick_scan(
    e3_structure="e3.pdb",
    glue_structure="glue.sdf",
    target_structures=["target1.pdb", "target2.pdb"],
    similarity_threshold=0.4
)
```

## Performance Tips

1. **Parallel Processing**: Set `n_workers` for multi-threaded scanning
2. **Pre-filter**: Use high similarity threshold initially, then relax
3. **Selective Docking**: Only dock top candidates
4. **Database Caching**: Pre-compute proteome surfaces

## References

- Molecular glue mechanisms: Słabicki et al., Nature 2020
- CRBN degraders: Krönke et al., Science 2014
- Pharmacophore methods: Wolber & Langer, J Chem Inf Model 2005
- Protein surface analysis: Sael & Kihara, Proteins 2012

## License

Part of the molecularGlues project.
