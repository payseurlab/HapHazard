# HapHazard

The HapHazard genomic ancestry simulator 

## Using HapHazard 

**For complete installation instructions, please see the manual (`Users_Manual.txt` in this repository).**

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

**In main folder: These are the files needed to compile and run HapHazard (see `Users_Manual.txt`). All other folders are unecessary for using HapHazard.** 

In folder `data_processing_scripts`, the following scripts can be found: 
- `Sim_Extract.R` : Will extract all simulation files from a large tar file of output from simulations run in parallel. 
- `Sim_Summary.R` : Will summarize junction density distributions from across the simulation output.
- `HH_DataProcDriver.R` : Can be used to coordinate processing of simulation output (run in parallel) using the scripts below.
- `HH_analysis_functions.R` : Contains functions used by the scripts below.
- `HH_ExtractExperiment.R` : Will extract all simulation files from a large tar file of output from simulations run in parallel. 
- `HH_BlockLength.R` : Will create a dataframe containing block length measurements. 
- `HH_JunctionDensity.R` : Will create a dataframe containing junction densities in 1 cM windows. 

In folder `run_sims_on_Condor`, the following scripts can be found: 
- `DAG-o-RAMA.pl` : creates a file structure to submit a number of parallel jobs on Condor. 
- `DAG_jobs.txt` : file that designates the jobs you would like to be run in parallel. 
- `MakeSimInp.pl` : file used by DAG-o-RAMA.pl to make the input file for the program (analogous to MakeHapHazInp.pl in the main folder)
- `HH_executable.sh` : Condor executable file
- `shared` : folder containing additonal optional processing scripts that can be run with the simulator.

Files in the folder `data_processing_scripts` can be used for post-processing of simulations, but are tailored to the file setup specific to the way we ran the simulator using Condor. Files in the folder `run_sims_on_Condor` are used to set up such a run. These files are used by putting them all in a folder with the HH_template.inp, modifying `DAG_jobs.txt` to match the simulations you wish to run, and then running `perl DAG-o-RAMA.pl`. This will generate a Condor-format DAG and a file structure. **Again, these scripts will not be applicable for other systems, but are provided as a potential basis for customization to other systems.** 
