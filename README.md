# Kupffer cells are key mediators of antigen-specific immune tolerance by peptide-coupled red blood cells
This repository contains the code to reproduce the scRNA-seq results from the paper _Kupffer cells are key mediators of antigen-specific immune tolerance by peptide-coupled red blood cells_ by Kalaitzaki et al. (currently under review). 

To generate the paper figures:
* Download the pre-processed Seurat objects from [Zenodo](https://doi.org/10.5281/zenodo.22029892) and place the `.rds` in the `processed_data` folder.
* Run the `generate_figures_all_cells.R` and `generate_figures_myeloid.R` scripts.

If you wish to do so, the code in `data_preprocessing.R` re-generates the pre-processed data and expects the Cellbender output matrices to be downloaded in the `raw_data` folder. These can be found on GEO (accession number [GSE342151](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE342151)).
