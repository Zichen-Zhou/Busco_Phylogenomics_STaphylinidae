# BUSCO Phylogenomics of Staphyliniformia 

Author and maintainer: Zichen Zhou.

Code accompanying the manuscript:

> *Reconstructing a highly supported backbone tree for the rove beetles (Staphylinidae) using BUSCO loci.*  
> Zichen Zhou & Alfried P. Vogler. Molecular Phylogenetics and Evolution (in review, 2026).

This repository contains the analysis code used for a BUSCO-based phylogenomic study of Staphyliniformia, including sequence assembly, BUSCO locus extraction and alignment, non-ortholog detection, locus evaluation, supermatrix construction, phylogenetic inference, topology tests, and node-support/concordance analyses.

---

## Pipeline overview

```text
01_GenomeAssembly
        ↓
02_BUSCOextraction
        ↓
03_LocusEvaluation
        ↓
04_Supermatrix
        ↓
05_PhylogeneticAnalysis
```

## Repository structure

```text
.
├── README.md
├── DEPENDENCIES.md
├── LICENSE
│
├── 01_GenomeAssembly/
│   ├── 1-script_nuclear_assembly_WGS.sh
│   ├── 1-script_nuclear_assembly_RNASeq.sh
│   └── 1-script_nuclear_assembly_EST.sh
│
├── 02_BUSCOextraction/
│   ├── 2.1-code_BUSCO_extraction.sh
│   └── 2.2-code_alignment.sh
│
├── 03_LocusEvaluation/
│   ├── 3.1-code_genetree_building.sh
│   ├── 3.2-code_non-ortholog_detection&removal.sh
│   ├── 3.2-script_apply_treeshrink_removal_4.0.py
│   ├── 3.3-script_alignment_for_locus_evaluation.sh
│   └── 3.4-code_locus_evaluation.sh
│
├── 04_Supermatrix/
│   ├── 4.1-code_supermatrix_loci_list.R
│   ├── 4.2-code_supermatrix_concatenation.sh
│   └── 4.2-script_filter_by_completeness.sh
│
└── 05_PhylogeneticAnalysis/
    ├── 5.1-code_phylogenetic_analysis.sh
    ├── 5.1-script_extract_main_GHOSTtree.py
    ├── 5.2-code_topology_test_RFdistance.R
    ├── 5.3-code_astral_tree_concensus.R
    ├── 5.4-code_topology_test_AUtest.sh
    └── 5.5-code_node_support_comparison.R
```

## Pipeline-to-paper map

| Analysis step | Code | Manuscript section |
|---|---|---|
| Genome and transcriptome assembly | `01_GenomeAssembly/` | Sequencing, assembly, and BUSCO recovery |
| BUSCO locus extraction and initial alignment | `02_BUSCOextraction/` | Sequencing, assembly, BUSCO recovery, and alignment |
| Non-ortholog removal and locus evaluation | `03_LocusEvaluation/` | Non-ortholog cleaning and locus evaluation |
| Locus selection and supermatrix construction | `04_Supermatrix/` | Filtering strategies and dataset construction |
| Concatenation- and coalescent-based phylogenetic analyses | `05_PhylogeneticAnalysis/5.1-*` | Phylogenetic reconstruction |
| RF-distance comparisons and ASTRAL consensus analyses | `05_PhylogeneticAnalysis/5.2-*`, `5.3-*` | Phylogenetic reconstruction |
| AU tests and constrained analyses | `05_PhylogeneticAnalysis/5.4-*` | Topology testing and constrained analyses |
| Node-support and concordance comparisons | `05_PhylogeneticAnalysis/5.5-*` | Concordance analysis |

## Data

The datasets associated with this study are archived on Zenodo:

- Data archive: `[DATA ZENODO DOI TO BE ADDED]`
- Archived release of this code repository: `[CODE ZENODO DOI TO BE ADDED]`

## Software

Software and package dependencies are listed in [`DEPENDENCIES.md`](DEPENDENCIES.md).

## Usage

Files labelled `code` are collections of commands used during the corresponding analysis step. They are retained as records of the analysis workflow and are intended to be run section by section where appropriate.

Files labelled `script` are standalone helper scripts used by the corresponding analysis steps.

Personal working paths, environment activation commands, and cluster-specific settings not required for the analyses were removed from the public version. Analytical commands, model settings, filtering thresholds, and the order of the main processing steps were retained.

The internal dataset labels `taxa50`, `taxa75`, and `taxa90` are retained in several filenames because they were used during the analyses. These correspond to `matrix50`, `matrix75`, and `matrix90` in the manuscript.

## Citation

Please cite the associated manuscript when using this code:

```bibtex
@article{Zhou2026Staphyliniformia,
  title   = {Reconstructing a highly supported backbone tree for the rove beetles (Staphylinidae) using BUSCO loci},
  author  = {Zhou, Zichen and Vogler, Alfried P.},
  journal = {Molecular Phylogenetics and Evolution},
  year    = {2026},
  doi     = {TO BE ADDED}
}
```

## License

MIT — see [`LICENSE`](LICENSE).
