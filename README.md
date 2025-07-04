# Multi-Ancestry GWAS of Cardiovascular Traits

This repository contains the analysis scripts and code used for the manuscript: **"Multi-trait GWAS Reveals Divergent Local and Polygenic Architecture of Cardiovascular Traits Across African and European Ancestries in the UK Biobank"**.

## 🧬 About The Project

This study performs a comprehensive genetic analysis of four cardiovascular traits—systolic blood pressure (SBP), diastolic blood pressure (DBP), pulse rate, and maximum heart rate (MHR)—in 459,327 European (EUR) and 6,654 African (AFR) ancestry individuals from the UK Biobank.The goal is to characterize and compare the genetic architecture of these traits across diverse ancestries, highlighting both shared and population-specific genetic determinants.

Our analysis pipeline integrates multi-trait GWAS, local genetic correlation analysis, functional annotation, and novelty assessment to provide a deep dive into the genetic landscape of cardiovascular health.

### Key Findings
* **Novel Variant Discovery:** Identified 957 novel variants in the EUR cohort and 45 novel variants in the AFR cohort.
* **Divergent Architecture:** Revealed a profound divergence in the pleiotropic architecture of blood pressure. We identified 181 genomic loci with significant local genetic correlation between SBP and DBP in the European sample, whereas such signals were completely absent in the African cohort.
* **Allelic Heterogeneity:** Found high replication of genomic *risk loci* between ancestries but low replication of the specific *lead SNPs*, suggesting different causal variants may exist within the same functional genomic regions.
* **Population-Specific Enrichment:** Uncovered distinct functional pathways, with AFR-associated genes enriched for epigenetic processes and EUR-associated genes enriched for cardiac development and function.

---

## 🛠️ Analysis Pipeline Overview

The project employs a multi-stage analysis pipeline using a combination of custom and publicly available software.

* **Phenotypic Analysis:** Comparison of cardiovascular and anthropometric traits between ancestry groups was performed in **MATLAB**.
* **Multi-Trait GWAS:** Joint analysis of the four cardiovascular traits was conducted using **JASS** to boost statistical power and identify pleiotropic loci.
* **Local Genetic Correlation:** Local SNP-heritability and genetic correlations were estimated using **LAVA** to investigate shared genetics at the locus level.
* **Functional Annotation:** Gene mapping, functional consequence prediction, and enrichment analysis were performed using **FUMA** and **Enrichr**.
* **Novelty Assessment:** Lead SNPs were checked against the GWAS Catalog using **LDLink** to distinguish between novel and known associations.
* **Heritability Estimation:** Genome-wide heritability and genetic correlations were estimated using **LDSC**.

---

## 💾 Data Availability

The primary individual-level phenotype and genotype data used in this study were obtained from the **UK Biobank Resource**. These data are subject to a Material Transfer Agreement and are not publicly available. Access can be obtained by approved researchers through a formal application to the UK Biobank.

* **UK Biobank Access Application:** [https://www.ukbiobank.ac.uk/enable-your-research/apply-for-access](https://www.ukbiobank.ac.uk/enable-your-research/apply-for-access)

Publicly available reference data used in our analyses can be downloaded from the following sources:
* **1000 Genomes Project:** [https://www.internationalgenome.org/data](https://www.internationalgenome.org/data)
* **GWAS Catalog:** [https://www.ebi.ac.uk/gwas/](https://www.ebi.ac.uk/gwas/)
* **LDSC LD Scores:** [https://alkesgroup.broadinstitute.org/LDSCORE/](https://alkesgroup.broadinstitute.org/LDSCORE/)
* **LAVA Genomic Loci:** [https://ctg.cncr.nl/software/lava/](https://ctg.cncr.nl/software/lava/)

---

## ✍️ How to Cite

If you use the code or findings from this project, please cite our manuscript:

> Sinkala, M., Elsheikh, S., Mbiyavanga, M., Cullinan, J., & Mulder, N. (2024). Multi-trait GWAS Reveals Divergent Local and Polygenic Architecture of Cardiovascular Traits Across African and European Ancestries in the UK Biobank.

## 💻 Repository Contents

This repository contains the core scripts for reproducing the analysis.

* `final_HPC_heartpipeline.m`: The main MATLAB script for phenotypic analysis, data integration, and generating figures.
* `create_jass_input.m`: MATLAB script to format GWAS summary statistics for JASS.
* `get_LDlink_GWAS_Catalog.m`: MATLAB script to perform novelty assessment using the LDLink API.
* `run_ldsc_local.sh`: BASH script to run LD Score Regression for heritability and genetic correlation.
* `parallel_cluster_lava_analysis.R`: R script to run LAVA for local genetic correlation analysis on an HPC cluster.
* `process_mtag_to_fuma.sh`: BASH script for pre-processing multi-trait analysis output.
* `run_jass_project.txt`: Example command for executing the JASS multi-trait analysis.
