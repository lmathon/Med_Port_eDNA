# Benchmarking biodiversity of seaports with eDNA and nearby marine reserves

This repository contains the codes to analyze eDNA data collected in ports and open areas in the Mediterranean French coast.


### Metadata
 
The [Metadata](00_Metadata/) folder contains the field metadata for each sample : sample code, date, station, site, project, sample type, volume, depth of sampling, habitat and coordinates.


### Analyses on Teleo barcode

The [01_Analyses_teleo](01_Analyses_teleo/) folder contains all the scripts to clean the data and run the ecological analyses.

#### Data

The [00_data](01_Analyses_teleo/00_data/) folder contains all the data files, matrices of sequences per sample, for the differnet sampling campaign.

The [01_Format_data](01_Analyses_teleo/01_Format_data/) folder contains the scripts to clean and assemble this data.


#### Analyses


The [02_Analyses](01_Analyses_teleo/02_Analyses/) folder contains the scripts to run the ecological analyses on the teleo eDNA data.
The folder [Analyses_reduced_species_list](01_Analyses_teleo/02_Analyses/Analyses_reduced_species_list/) contains teh same analyses scripts but for the matrices containing only species validated by experts.


#### Outputs and plots

The [03_Outputs](01_Analyses_teleo/03_Outputs/) folder contains the figures and tables from the analyses scripts.
The [04_Plots](01_Analyses_teleo/04_Plots/) folder contains the scripts to build the final figures for the paper, and the figures and tables.


