## Code for "Systematic in vitro evolution in Plasmodium falciparum reveals key determinants of drug resistance"

This repository contains original code used to analyze data in "Systematic in vitro evolution in Plasmodium falciparum reveals key determinants of drug resistance" by Luth et al.

- CNV_calling_validation contains shell scripts/Python scripts and metadata (lists of sample .bam filenames) used to call and validate CNVs. DELLY (https://github.com/dellytools/delly) calls were used in the validation step.
- dNdS_analysis contains Python scripts that calculate "expected" nonsnonymous/synonymous SNVs per gene, used to approximate dN/dS in the study's dataset of in vitro evolved mutations as well as the Pf6 dataset of field variants
- Under Jupyter_notebooks:
    - final_resistome_CNVs.ipynb contains code used to create Figure 2
    - SNV_evolved_vs_Pf6_dNdS.ipynb contains code used to create Figure 3D
    - CNV_visualizer.ipynb contains analyses that informed the heuristic CNV validation method
- reference_files contains the P. falciparum 3D7 reference genome used in this study, list of gene intervals, and PON (panels of normals) used to denoise read counts for CNV calling

