# E. coli Pangenome-Scale Metabolic Model Reconstruction and Analysis: "Pangenome-Scale Metabolic Network Reconstruction Reveals a Diverse Genetic Basis for Metabolic Reactions"

This repository contains data and code for the pangenome-scale metabolic reconstruction of *Escherichia coli*, enabling detailed analysis of its metabolic diversity. By analyzing the gene-to-protein-to-reaction associations (GPRs) across thousands of strains, this study explores the genetic basis of E. coli's metabolism and reveals the evolutionary dynamics shaping its metabolic functions.



## Method Overview (Figures)

<table>
  <tr>
    <td width="50%">
      <img src="docs/workflow.png" alt="Automated pipeline overview" width="100%">
      <p align="center"><em>rare essential genes of <i>E. coli</i>: a BiGG data analysis.</em></p>
    </td>
    <td width="50%">
      <img src="docs/panel1.jpg" alt="Pangenome annotation and rare genes analysis" width="100%">
      <p align="center"><em>pangenome annotation and rare genes analysis.</em></p>
    </td>
  </tr>
</table>





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

