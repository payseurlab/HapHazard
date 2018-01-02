# HapHazard
The HapHazard genomic ancestry simulator

Dependencies:
GNU scientific library: last compiled with v2.2 downloaded from http://ftp.wayne.edu/gnu/gsl/
Perl 5 (or later version)

Installation:
1. Download the zip file from this repsoitory.
2. Navigate to the folder containing the zip in the terminal and unzip the file using the command:
    $ unzip HapHazard-master.zip
3. Compile the program from source using the g++ compiler and linking GSL using the command
    $ g++ -L/usr/local/lib main.cpp -o HapHazard -lgsl -lgslcblas -lm
    The link command -L/usr/local/lib links the compilation to gsl's default installation. If you've chosen to install it
    elsewhere change this command appropriately. This command should produce an executable file named 'HapHazard'

Setting your First (Test) Run

The installation comes with the file 'HH_template.inp' which serves as an input template. The default parameters in it
can be used to run a test simulation.

1. Install Perl of you have not already, and run the input generator script by typing:
   $ perl MakeHapHazInp.pl HH_template.inp
   This command will run the input generation by processing the parameters listed in "HH_templatr.inp" and formatting them into a list that HapHazard can process which is stored in a file called "Homer.inp" The file is named "Homer.inp" because "Homer" is the experiment name parameter given in "HH_template.inp". Variables and default parameters for HapHazard are listed in "HH_template.inp" with brief explanations and instructions for how to use them.
   
2. To run Haphazard (on Linux or Unix) type:
   $ ./HapHazard Homer.inp 0 555
   On windows, simply omit the "/." symbol. HapHazard requires command line parameters to run. The first of the name of the input file as the first command line parameter, here it's "Homer.inp". This will vary depending on how you name your experiments. Next, you must specify a simulation number, in this case "0", and a random seed number, here it's "555". You can change these last two, but make sure that no two simulations have the same number or random seed. Having the same number will cause previous results to be overwritten, and any two simulations with exactly the same seed numbers (and parameter values) will be exact copies of each other. 
   To verify that the simulation is running properly check the output. The first line should say "HapHazard v1.0". The output should look like the contents of the "HapHazard_test_output.txt" file in this repository.
   
3.After finishing, their should be a folder called "Homer_0" in the HapHazard-master directory. In it, there should be the files. "Homer.desc" list the parameters and their values that were used in the experiment. "Homer.0.499.clines" contains a list of the counts of each ancestry by chromosome and marker position, for each deme. The name of the file indicates that these were experiment Homer, simulation number 0, sample from generation 499, and the data are for measuring clines indicated by the extension. The file ending with ".gclines", contains the same but the congress are genotypes. There should also be three folders, "BL_CHR0", "BL_CHR1", and "BL_CHR2".
