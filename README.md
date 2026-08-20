# Kupffer cells are key mediators of antigen-specific immune tolerance by peptide-coupled red blood cells
This repository contains the code to reproduce the scRNA-seq results from the paper _Kupffer cells are key mediators of antigen-specific immune tolerance by peptide-coupled red blood cells_ by Kalaitzaki et al. (currently under review). 

The scripts to generate the paper figures require the pre-processed data which can be downloaded from Zenodo: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22029892.svg)](https://doi.org/10.5281/zenodo.22029892). The `.rds` files should be placed in the `processed_data` folder.

Alternatively, the code in `data_preprocessing.R` re-generates the processed data and expects the Cellbender output matrices to be downloaded in the `raw_data` folder. These can be found on GEO (accession number [GSE342151](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE342151)).
