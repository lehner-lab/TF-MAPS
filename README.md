# TF-MAPS
TF-MAPS: fast high-resolution functional and allosteric mapping of DNA-binding proteins. 
Here you'll find source code to reproduce the computational analyses in the [bioRxiv paper](https://www.biorxiv.org/content/10.1101/2025.10.20.683418v1). 

## Required Data
The following data as RData formats are available for download from [here zeondo](https://zenodo.org/records/17457956), and copied to your "base_dir" folder for running the analysis.    
(1) the read counts (DiMSum outputs for DNA-binding assay "dimsum_output_b1h.RData" & "dimsum_output_ss.RData"),   
(2) the combined score data with all annotations for all amino acid substitutions ("combined_muts_hn.RData", "combined_muts_fg.RData", "combined_muts_fp.RData"),   
(3) the positional aggregated data ("positional_hn.RData", "positional_fg.RData",  "positional_fp.RData"),  
(4) Combined data for FOXG1 and FOXP1 for direct comparisons ("muts_fg_fp_comparisons.RData", "positional_fg_fp_comparisons.RData" ) 

## System Requirements
R (GGally, ggplot2, ggpubr, gplots, stringr,dplyr) 

## Information about the provided R scripts 
1.	fitting_loess.R  
This function fits the relationship between DNA-binding assay measurements and protein abundance for variants. It allows calculation of binding residuals, representing variant effects on DNA binding independent of protein abundance.  
2.	distance_decay.R  
This function computes the relationship between the absolute binding residuals of variants and their distance to DNA, as determined from PDB structures. The resulting distance-corrected binding residuals quantify the anisotropy of DNA-binding effects.  
3.	chimera_defattr_visualisation.R  
This script generates DEFATTR files for ChimeraX, enabling visualisation of median variant effects on the 3D structure of the protein.  
4.	abundance_ic50_mean_abundance_comparisons.R  
This script calculates IC50 values from spectinomycin assays using DiMSum output across 21 concentrations. These IC50 values are compared with mean enrichment scores from three selected concentrations to guide the choice of scoring matrices.

Usage: Example commands and instructions for each function are provided within the respective R scripts alongside the function definitions.  

## Additional Information 
To reproduce the part from Illumina sequencing reads to DiMSum outputs, please use DiMSum v1.3.2. Download the FastQ files from European Nucleotide Archive (ENA) with accession number [PRJEB97482](https://www.ebi.ac.uk/ena/browser/view/PRJEB97482) to your base directory. Parameters for the DiMSum run are provided in the [manuscript] (https://www.biorxiv.org/content/10.1101/2025.10.20.683418v1).   
