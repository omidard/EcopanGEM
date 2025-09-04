# E. coli Pangenome-Scale Metabolic Model Reconstruction and Analysis: "Pangenome-Scale Metabolic Network Reconstruction Reveals a Diverse Genetic Basis for Metabolic Reactions"

This repository contains data and code for the pangenome-scale metabolic reconstruction of *Escherichia coli*, enabling detailed analysis of its metabolic diversity. By analyzing the gene-to-protein-to-reaction associations (GPRs) across thousands of strains, this study explores the genetic basis of E. coli's metabolism and reveals the evolutionary dynamics shaping its metabolic functions.



## Method Overview (Figure)

<p align="center">
  <img src="docs/workflow.png" width="900" alt="Automated pipeline: genomes → QA/QC → annotation (BAKTA) → pangenome (CD-HIT) → GEM reconstruction (CarveMe) → panGPRs → neighborhood analysis; with protein stoichiometry, 3D modeling, and structural analysis integration.">
</p>

*Figure.* rare essential genes of E.coli: A bigg data analysis.

<p align="center">
  <img src="docs/panel1.jpg" width="900" alt="Automated pipeline: genomes → QA/QC → annotation (BAKTA) → pangenome (CD-HIT) → GEM reconstruction (CarveMe) → panGPRs → neighborhood analysis; with protein stoichiometry, 3D modeling, and structural analysis integration.">
</p>

*Figure.* pangenome annotation and rare genes analysis.




## Overview

We constructed strain-specific metabolic models for 2,377 fully sequenced *E. coli* strains, covering approximately 2,700 metabolic reactions. The resulting pangenome-scale model, or "panGEM," connects metabolic genotypes to phenotypes, offering insights into the genetic diversity and evolutionary history of *E. coli* metabolism.


## Repository Structure

- **data/**  
  Contains data files in csv format.
  
- **scripts/**  
  Contains code for data processing, pangenome model reconstruction, and analysis.

- **notebook/**  
  Jupyter notebook for exploratory data analysis, visualization, and results replication.


Files required for running the notebook are deposited on zenodo: https://zenodo.org/records/14028473
GEMs are deposited on Zenodo:  https://zenodo.org/records/13825392

