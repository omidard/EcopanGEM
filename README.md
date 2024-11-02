# E. coli Pangenome-Scale Metabolic Model Reconstruction and Analysis: "Pangenome-Scale Metabolic Network Reconstruction Reveals a Diverse Genetic Basis for Metabolic Reactions"

This repository contains data and code for the pangenome-scale metabolic reconstruction of *Escherichia coli*, enabling detailed analysis of its metabolic diversity. By analyzing the gene-to-protein-to-reaction associations (GPRs) across thousands of strains, this study explores the genetic basis of E. coli's metabolism and reveals the evolutionary dynamics shaping its metabolic functions.

## Overview

We constructed strain-specific metabolic models for 2,377 fully sequenced *E. coli* strains, covering approximately 2,700 metabolic reactions. The resulting pangenome-scale model, or "panGEM," connects metabolic genotypes to phenotypes, offering insights into the genetic diversity and evolutionary history of *E. coli* metabolism.

### Key Findings

1. **Genetic Basis of Metabolic Properties**  
   The GPRs derived from these models reveal a species-wide genetic basis for specific metabolic properties, showcasing genetic diversity and conservation in E. coli metabolism.

2. **Diversity in Core Metabolism**  
   Many rare genes associated with core metabolic pathways come from diverse protein families, demonstrating both gene fragmentation and horizontal gene transfer.

3. **Genomic Neighborhood Variability**  
   Rare genes show high variability in genomic neighborhoods, often associated with transposable elements, suggesting a dynamic evolutionary process.

4. **Amino Acid Pathway Enrichment**  
   Pathways for aromatic amino acids and branched-chain amino acids are enriched with rare genes, highlighting complex evolutionary pressures.

## Repository Structure

- **data/**  
  Contains data files in csv format.
  
- **scripts/**  
  Contains code for data processing, pangenome model reconstruction, and analysis.

- **notebook/**  
  Jupyter notebook for exploratory data analysis, visualization, and results replication.




