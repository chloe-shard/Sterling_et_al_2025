# Sterling_et_al_2025

Overview

Code supporting the publication: Mechanistic deconvolution reveals DHODH as the key target of the KDM4 inhibitor QC6352, Nature Chemical Biology, 2026.
Author: Jayden R. Sterling, Jacinda R. Holtsmark, George L. Joun , Karen Tran, Tian Y. Du, Mani Kuchibhotla, Terrance G. Johns, Antoine deWeck, Lipin Loo, G. Gregory Neely, Chloe Shard, Guillermo Gomez, Ashwini Patil, Diana D. Shi, Julie-Aurore Losman, Angus Saxton, Jennifer R. Baker, Zhihe Lei, Matthew Holland, Paul E. Brennan, Mathew Graus, Mathias Francois, Andreas Kraemer, Susanne Mueller, Yuchen Feng, Paul Workman, Lenka Munoz

- Chronos_CRISPR_analysis.ipynb uses machine learning to train a model to identify CRISPR gene effects over time.
- Limma_CRIPSR_analysis.ipynb uses linear models and empirical Bayes statistical analysis as a supplementary method for identifying CRISPR gene effects over time. 
- RNA-seq_volcano_plots.ipynb creates volcano plots for RNA-sequencing data.
- GO analysis.ipynb analyses and visualises RNA-sequencing differentially expressed genes.
- Bulk_seq_GSC_KDM_heatmap_Rscript creates a heatmap of KDM family expression across GSC from GSE119834 (PMID:30948495).
- scRNAseq_GBM_malignant_cells_KDM_heatmap_Rscript creates a data matrix of KDM family expression in malignant cells from scRNAseq data. the output data matrix should then be imported into https://software.broadinstitute.org/morpheus/ to create a heatmap. scRNAseq is from patient biopsies from 6 public scRNAseq datasets: GSE131928, 20 patients, GSE173280, 10 patients, GSE182109, 8 patients, GSE141383, 8 patients, EGAD00001006206, 13 patients and Ebert et al., 2020, 3 patients. scRNAseq data from Ebert et al., 2020 can be requested by contacting G.Gomez (guillermo.gomez@unisa.edu.au)

Requirements

All required packages are listed in the import cells of each notebook.

Python 3 environment: Chronos_CRISPR_analysis.ipynb and RNA-seq_volcano_plots.ipynb

R environment: all remaining notebooks


How to Run

Download or clone the repository.

Open the notebook you want to run.

Ensure you're using the correct environment (Python or R kernel) and all packages are up to date.

Run all cells.


Citation

If you use this code, please cite the associated publication.
