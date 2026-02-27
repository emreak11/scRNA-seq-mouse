**Exploratory preprocessing of E8.5 mouse endothelium single-cell RNA-seq data**

**Project Overview**
•This project is a small exploratory exercise focused on preprocessing single-cell RNA-seq data from E8.5 Mus musculus endothelium.
•The goal was not to answer a specific biological question, but to understand how raw count data is processed before downstream analysis.

**What I Did**
•Loaded raw single-cell RNA-seq count matrices
•Resolved gene identifier mapping using Ensembl annotation
•Identified mitochondrial genes
•Computed per-cell mitochondrial fraction (mt-frac)
•Applied quality-control filtering based on mt-frac
•Performed library-size normalization
•Applied log(1 + x) transformation

**Key Learning Points**
•Understood the importance of gene ID vs gene name mapping
•Learned how mitochondrial gene content is used for quality control
•Gained practical experience in preprocessing workflows for single-cell data
•Developed basic data validation and debugging skills during pipeline construction

**Technical Details**
•Implemented in Python
•Used pandas and basic numerical operations
•Manual threshold inspection for quality control

**Limitations**
•No PCA, UMAP, or clustering performed
•No biological interpretation step
•Focused only on preprocessing and data preparation
