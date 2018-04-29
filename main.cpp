/*
    *** HAP_HAZARD v1.0 ***

    ** NAMING CONVENTIONS **
    constants: ALL_CAPS
    classes: Capitalized Camel case. ex: ClassName
    methods: Camel case. ex: methodName()
    static variables: snake case beginning with "s_".  ex: s_static_variable
    private/protected variables: snake case, begins with "m_".  ex: m_private_or_protected_variable
    public variables: snake case, begins with "p_". ex: p_public_variable
    global variables: snake case, begins with a "g_". ex: g_global_variable
*/


// From the standard library
#include <iostream>
#include <vector>
#include <ctime>
#include <deque>

// from the GNU scientific library (GSL)
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// classes and methods used in this project
#include "junction.h"           // updated
#include "cnode.h"              // updated
#include "chromosome.h"         // updated
#include "autosome.h"           // updated
#include "gamete.h"             // updated
#include "sex_chromosome.h"     // updated
#include "cytoplasm.h"          // updated
#include "individual.h"         // updated
#include "gene.h"               // updated
#include "jungen_acc.cpp"       // updated
#include "landscape.h"          // updated
#include "interaction_graphs.h" // updated
#include "deme.h"               // updated
#include "metaPop.h"            // updated
#include "marker.h"             // updated
#include "experiment.h"         // updated

// CLASS STATIC VARIABLE DECLARATIONS
// class random number generators -- all will use the same RNG
// Comments for individual variables are in their class definitions
gsl_rng* Experiment::s_rand_num;
gsl_rng* MetaPop::s_rand_num;
gsl_rng* Deme::s_rand_num;
gsl_rng* Individual::s_rand_num;
gsl_rng* Chromosome::s_rand_num;

// Junction static variables
std::deque<Junction*> Junction::s_junction_pool;
int Junction::s_num_locations;
int Junction::s_num_recycled;
int Junction::s_num_junctions;

// Cnode static variables
int CNode::s_num_cnode;

// Chromosome static variables
double Chromosome::s_par_boundary;
double Chromosome::s_window_size;
int Chromosome::s_num_anc;

// Gamete static variables
std::deque<Gamete*> Gamete::s_gamete_reservoir;
int Gamete::s_num_gametes;

// Individual static variables
int Individual::s_num_chr;
int Individual::s_interference_opt;
int Individual::s_num_of_phenotypes;
std::vector<Gamete*> Individual::s_ind_gamete_pool;
std::vector<Landscape*> Individual::s_landscapes;
MatingSys Individual::s_gp_type;
double Individual::s_min_block_size;

// Deme static variables
MatingSys Deme::s_mat_sys;
int Deme::s_num_anc;

// Gene static variables
int Gene::s_num_anc, Gene::s_num_chr;

// Marker static variables
int Marker::s_num_demes;
int Marker::s_num_ancestries;
int Marker::s_mark_sample_size;
std::vector<std::string> Marker::s_genotype_keys;
std::ofstream Marker::s_cline_file;
std::ofstream Marker::s_g_cline_file;

// IntGraph static variables
int IntGraph::s_num_ancestries;

// Experiment static variables

// Object Recycling variables
std::deque<Junction*> Junction::s_junction_reservoir;
std::deque<CNode*> CNode::s_cnode_reservoir;
std::deque<Chromosome*> Chromosome::s_chromosome_reservoir;
std::deque<Individual*> Individual::s_individual_reservoir;

// Global variables
int g_generation;           // the current generation of the simulation, global because nearly every object/method uses it
std::ofstream g_check_file;   // used to check output for errors when output is too large to read in the terminal, needs to be called anywhere like std::cout

int main(int argc, char* argv[])
{
    // Just print the name to the terminal
    std::cout << "HAPHAZARD v1.0" << std::endl;

    // Open the check file and do the same
    g_check_file.open("check_file.txt");
    g_check_file << "HAPHAZARD v1.0" << std::endl;

    // Construct the experiment
    const char* INPUT_FILE_NAME = argv[1];       // the first arg is the input file's name
    const char* EXPERIMENT_NUMBER = argv[2];     // the second is the number for the experiment
    const char* RNG_SEED = argv[3];              // the third is the seed value for the random number generator
    Experiment hap_hazard(INPUT_FILE_NAME, EXPERIMENT_NUMBER, RNG_SEED);   // pass the args to the experiment constuctor

    // Run the experiment
    hap_hazard.experimentRun();

    // After the experiment is finished use the block data to make marker files -- requires the Perl script 'happhazard_makemarkers.pl'
    std::cout << "Making Marker Files..." << std::endl;
    system("perl haphazard_makemarkers.pl");

    // generate some R scripts to hold data for future analyses
    // we'll make a stringstream and pass strings formatted as R code to it
    // then print it to a file
    std::stringstream r_com;
    std::stringstream r_vars;

    // add the experiment name
    std::string rvar_exp_name = hap_hazard.getExperimentName();
    r_vars << "expName = \"" << rvar_exp_name << "\"" << std::endl;
    r_vars << "jobID = \"" << hap_hazard.getJobID() << "\"" << std::endl;

    // add the chromosome lengths
    r_vars << "chrLength=c(";
    std::vector<double> rvar_chromosome_lengths = hap_hazard.getChrLengths();
    std::vector<double>::iterator iter_cl;
    for(iter_cl = rvar_chromosome_lengths.begin() ; iter_cl < rvar_chromosome_lengths.end() - 1 ; iter_cl++ )
    {
        r_vars << (*iter_cl) << ",";
    }
    r_vars << (*iter_cl) << ") " << std::endl;

    // add the numbers of chromosomes -- i.e. indices in an array of chromosomes
    int rvar_num_chr = hap_hazard.getNumChr();
    r_vars << "chromosomes=c(" << std::endl;

    for( int c = 0 ; c < rvar_num_chr - 1 ; c++)
    {
        if( c == rvar_num_chr - 2 )     // because we're ignoring the cytoplasmic marker in the analyses that use this file
        {
            r_vars << c << ") " << std::endl;
        }
        else
        {
            r_vars << c << ",";
        }
    }

    // now we get the indicies for the list of summary generations
    std::vector<int> rvar_summary_generations = hap_hazard.getSummaryGenerations();
    r_vars << "generations=c(";
    std::vector<int>::iterator iter_sg;
    for(iter_sg = rvar_summary_generations.begin() ; iter_sg < rvar_summary_generations.end() - 2 ; iter_sg++ )
    {
        r_vars << (*iter_sg) << ",";
    }

    r_vars << (*iter_sg) << ") " << std::endl;

    // now the indices for the demes
    int rvar_number_of_demes = hap_hazard.getNumDemes();
    r_vars << "demes=c(";
    
    for( int d = 0 ; d < rvar_number_of_demes ; d++ )
    {
        if( d == rvar_number_of_demes )
        {
            r_vars << d << ") " << std::endl;
        }
        else
        {
            r_vars << d << ",";
        }
    }

    // the "junction demes" are the demes we are going to collect data for
    r_vars << "junction_demes=c(";
    
    for( int d = 0 ; d < rvar_number_of_demes ; d++ )
    {
        if( d == rvar_number_of_demes - 1 )
        {
            r_vars << d << ") " << std::endl;
        }
    }

    // get the window size/marker spacing
    r_vars << "winSize=" << hap_hazard.getWinSize() << std::endl;

    // add the stringstream to the file
    std::ofstream rvar;
    rvar.open("Rvars.R");
    rvar << r_vars.str();
    rvar.close();

    // close the check file
    g_check_file.close();

    return 0;
}
