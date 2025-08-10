# Morphological Analysis of Middle Pleistocene Skulls

This repository contains the code, data, and most results for the analyses presented in:

**Mishol N., Herzlinger G., Rak Y., Smilansky U., Carmel L., Gokhman D. (2025)**  
*Candidate Denisovan fossils identified through gene regulatory phenotyping*

## Overview

This repository includes:
- Preprocessed phenotype data and measurement tables
- Code used to generate all figures and statistical results
- Supplementary results and validation analyses
- Scripts to reproduce and extend the core analyses

## Reproducing the Main Analysis

1. Download the raw meaurements from Ni et al. 2021 (Morphobank)
2. process TNT file using proccess_TNT_file.py
3. Run add_metadata_and_change_synonyms.R
4. Run `Main_analyses.Rmd` in RStudio.
5. Outputs will be saved in the `results/` directory.