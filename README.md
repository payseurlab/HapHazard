# HapHazard
The HapHazard genomic ancestry simulator

Dependencies:
GNU scientific library: last compiled with v2.2 downloaded from http://ftp.wayne.edu/gnu/gsl/

Installation:
1. Download the zip file from this repsoitory.
2. Navigate to the folder containing the zip in the terminal and unzip the file using the command:
    $ unzip HapHazard-master.zip
3. Compile the program from source using the g++ compiler and linking GSL using the command
    $ g++ -L/usr/local/lib main.cpp -o HapHazard -lgsl -lgslcblas -lm
    The link command -L/usr/local/lib links the compilation to gsl's default installation. If you've chosed to install it
    elsewhere change this command appropriately. This command should produce and executable file named 'HapHazard'

Setting your First (Test) Run

The installation comes with the file 'HH_template.inp' which serves as an input template. The default parameters in it
can be used to run a test simulation.
