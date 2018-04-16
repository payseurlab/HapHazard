#!/bin/bash
tar -xzf R.tar.gz
export PATH=$(pwd)/R/bin:$PATH
tar -xzf SLIBS.tar.gz 
export LD_LIBRARY_PATH=$(pwd)/SS:$LD_LIBRARY_PATH
R CMD BATCH run_hap_hazard.R

