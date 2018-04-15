# HapHazard

The HapHazard genomic ancestry simulator 

## Using HapHazard 

**For complete installation instructions, please see the manual (Users_Manual.txt in this repository).**

System Requirements: 
- Linux (we used Ubuntu 14.04 and later) 
- Cygwin (if using Windows)
- 2GB of RAM (at least)

Dependencies: 
- GNU Scientific Library
- Perl 5 (or later version) 

Installation: 
1. Download the repository. 
2. Navigate to the folder in the terminal. 
3. Compile the program from source using the g++ compiler and linking GSL using the command: 
  `$ g++ -L/usr/local/lib main.cpp -o HapHazard -lgsl -lgslcblas -lm` 
 

## Contents 

This folder contains the files needed to compile HapHazard, as well as some extra files for processing. The additional scripts are specific to simulations run in parallel on HT Condor and will only be useful for similar systems or as a basis for customization to another system. 

In folder `data_processing_scripts`, the following scripts can be found: 
- `Sim_Extract.R` : Will extract all simulation files from a large tar file run in parallel. 
- `Sim_Summary.R` : Will summarize junction density distributions from across the simulation output. 

In folder `run_sims_on_Condor`, the following scripts can be found: 
- `DAG-o-RAMA.pl` : 
- `DAG_jobs.txt` : 
