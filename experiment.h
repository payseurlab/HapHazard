#ifndef EXPERIMENT_H
#define EXPERIMENT_H

#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>

#include <sstream>
#include <string>
#include <stdio.h>
#include <ctime>
#include <sys/stat.h>
#include <stdlib.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "junction.h"
#include "cnode.h"
#include "chromosome.h"
#include "autosome.h"
#include "sex_chromosome.h"
#include "cytoplasm.h"
#include "individual.h"
#include "gene.h"
#include "jungen_acc.cpp"
#include "landscape.h"
#include "interaction_graphs.h"
#include "deme.h"
#include "metaPop.h"
#include "marker.h"


class Experiment
{
    public:
        // Experiment constructor prototype
        Experiment(const char* iF, const char* jI, const char* rS);

        // Input Functions
        void experimentInput();
        void fileInput();

        // Getters
        std::string getExperimentName() { return m_experiment_name; }
        int getNumChr() { return m_number_of_chromosomes; }
        std::vector<double> getChrLengths() { return m_chromosome_lengths; }
        std::vector<int> getSummaryGenerations() { return m_summary_generations; }
        int getSimNum() { return m_sim_number; }
        const char* getJobID() { return m_job_id; }
        int getNumDemes() { return m_number_of_demes; }
        double getWinSize() { return m_marker_spacing; }
        int getSampSize() { return m_sample_size; }

        // Experiment Instantiation and Run-time
        void experimentSetup();
        void experimentRun();

        void genotypeSample(Individual**);
        void textMarkers(std::string, std::string);
        const char* makeFileName(std::string);
        const char* makeMarkerFileName(std::string);
        std::string makeBlockFileName(std::string, int);
        void makeMarkers(double s);

        void summary(Individual**);
        void dbSummary(Individual**);
        void blockSummary(Individual**);

        void makeDescFile();

        void memMonitor();
        void memClean();

        static gsl_rng * s_rand_num;

    private:
        int m_number_of_simulations;
        int m_sim_number;

	// Need Filenames for I/O
        // Input
        const char* m_INPUTFILE;

        const char* m_job_id;
        int m_number_of_generations;

        std::string m_experiment_name;
        //std::string m_job_id;
        int m_seed;

        int m_number_of_demes;
        std::vector<int> m_deme_sizes;
        std::vector<double> m_migration_rates;
        MatingSys m_mating_system;

        int m_number_of_ancestries;
        std::vector<double> m_ancestry_frequencies;

        int m_number_of_chromosomes;
        int m_interference_option;
        std::vector<double> m_chromosome_lengths;
        std::vector<ChrType> m_chromosome_types;
        double m_par_length;

        int m_number_of_phenotypes;
        std::vector<Gene*> m_genes;
        std::vector<IntGraph*> m_network;
        std::vector<double> m_exp_phen_values;
        std::vector<Landscape*> m_phen_landscapes;

        // Summary Stats and other Output
        int m_sample_size;
        std::vector<Marker*> m_markers;
        double m_marker_spacing;
        double* m_junction_density_table;
        std::vector<int> m_summary_generations;

        // MetaPop Parameters
        std::vector<Gamete*> m_gamete_pool;

        // The metapop
        MetaPop* m_hap_hazard;

        // Output
        // SETUP FILE STREAMS FOR OUTPUT
        // ___.clines, ___.gclines, ___.pca, ___.jden, ___.jfrq, ___.jdir, ___.bla
        bool m_clines_on;
        std::ofstream m_os_clines;

        std::ofstream m_os_gclines;
        bool m_gclines_on;

        std::ofstream m_os_pca;
        bool m_pca_on;

        std::ofstream m_os_junction_density;
        bool m_junction_density_on;

        std::ofstream m_os_junction_freq;
        bool m_jun_freq_on;

        std::ofstream m_os_junction_dir;
        bool m_jun_dir_on;

        std::ofstream m_os_blocks;
        bool m_blocks_on;
};

//Experiment constructor
Experiment::Experiment(const char* iF, const char* jI, const char* rS):
    m_INPUTFILE(iF), m_job_id(jI)
{
    // GSL RNG setup
    s_rand_num = gsl_rng_alloc(gsl_rng_taus);
    m_seed = strtol(rS,0,10) + strtol(jI,0,10);
    gsl_rng_set(s_rand_num, m_seed);

    fileInput();

    std::stringstream dir_name;
    dir_name << m_experiment_name << "_" << jI << "/";
    std::string sdir = dir_name.str();
    const char* dir = sdir.c_str();
    mkdir(dir, 0775);

    if (m_blocks_on)
    {
        for( int c = 0 ; c < m_number_of_chromosomes ; c++ )
        {
            dir_name.str("");
            dir_name << m_experiment_name << "_" << m_job_id << "/BL_CHR" << c << "/";
            sdir = dir_name.str();
            dir = sdir.c_str();
            mkdir(dir, 0775);
        }
    }
    makeDescFile();

    experimentSetup();
}

// MAKE A DESCRIPTION FILE FOR THE EXPERIMENT
void Experiment::makeDescFile()
{
    std::ofstream description;

    std::string appendix = "desc";
    std::stringstream filename;
    filename <<  m_experiment_name << "_" << m_job_id << "/" << m_experiment_name << "." << appendix;
    std::string desc_file = filename.str();
    std::cout << desc_file << std::endl;

    description.open(desc_file.c_str(), std::ios_base::trunc);

    description << "Experiment Name: " << m_experiment_name << std::endl;
    description << "Seed: " << m_seed << std::endl;
    description << "Number of Simulations: " << m_number_of_simulations << std::endl << std::endl;

    description << "DEMOGRAPHIC PARAMETERS" << std::endl << std::endl;

    description << "Number of Generations: " << m_number_of_generations << std::endl;
    description << "Number of Demes: " << m_number_of_demes << std::endl;
    description << "Number of Ancestries " << m_number_of_ancestries << std::endl;
    description << "Mating System: " << m_mating_system << std::endl;
    description << "Deme\tSize\t";

    for( int i = 0 ; i < m_number_of_ancestries ; i++ )
    {
        description << "Anc_" << i << "\t";
    }
    description << std::endl;

    for( int i = 0 ; i < m_number_of_demes ; i++ )
    {
        description << i << "\t\t" << m_deme_sizes[i];
        for( int j = 0 ; j < m_number_of_ancestries ; j++ )
        {
            description << "\t\t" << m_ancestry_frequencies[m_number_of_ancestries * i + j];
        }
        description << std::endl;
    }

    description << "Migration Rates: ";
    std::vector<double>::iterator iter_d;
    for( iter_d = m_migration_rates.begin() ; iter_d < m_migration_rates.end() ; iter_d++ )
    {
        description << " " << (*iter_d);
    }
    description << std::endl << std::endl << "GENOMIC PARAMETERS" << std::endl << std::endl;

    description << "Number of Chromosomes: " << m_number_of_chromosomes << std::endl;
    description << "Recombination Model: " << m_interference_option << std::endl;
    description << "Pseudoautosomal Region: " << m_par_length << std::endl << std::endl;

    description << "Chr\tType\tLength" << std::endl;

    for( int i = 0 ; i < m_number_of_chromosomes ; i++ )
    {
        description << i << "\t\t" << m_chromosome_types[i] << "\t\t" << m_chromosome_lengths[i] << std::endl;
    }
    description << std::endl << std::endl << "GENETIC PARAMETERS" << std::endl << std::endl;

    description << "Number of Phenotypes: " << m_number_of_phenotypes << std::endl;
    description << "Phenotype\tExp_Value" << std::endl;
    int i = 0;
    for(iter_d = m_exp_phen_values.begin() ; iter_d < m_exp_phen_values.end() ; iter_d++ )
    {
        description << i << "\t\t\t" << (*iter_d) << std::endl;
        i++;
    }

    description << std::endl << "Genes:" << std::endl << "Chr\tPos\tPhen";

    for(i = 0 ; i < m_number_of_ancestries ; i++ )
    {
        description << "\tAddEff_" << i ;
    }
    description << std::endl;

    std::vector<Gene*>::iterator iter_g;
    for(iter_g = m_genes.begin() ; iter_g < m_genes.end() ; iter_g++ )
    {
        description << (**iter_g).getChr() << "\t" << (**iter_g).getPos() << "\t" << (**iter_g).getPhenotype();
        std::vector<double> aff = (**iter_g).getAddFX();
        std::vector<double>::iterator iterF;

        for( iterF = aff.begin() ; iterF < aff.end() ; iterF++ )
        {
            description << "\t\t" << (*iterF);
        }
        description << std::endl;
    }

    description << "Interactions: " << std::endl;
    description << "G1_chr\tG1_pos\tAncOne\tG2_chr\tG2_pos\tAncTwo\tMaxEff\tPhen\tModel" << std::endl;
    std::vector<IntGraph*>::iterator iter_ig;
    for(iter_ig = m_network.begin() ; iter_ig < m_network.end() ; iter_ig++ )
    {
        description << (*((**iter_ig).getLocusA())).getChr() << "\t" << (*((**iter_ig).getLocusA())).getPos() << "\t" << (**iter_ig).getAncA()
        << "\t" << (*((**iter_ig).getLocusB())).getChr() << "\t" << (*((**iter_ig).getLocusB())).getPos() << "\t" << (**iter_ig).getAncB() << "\t"
        << (**iter_ig).getMaxSel() << "\t\t" << (**iter_ig).getPhenotype() << "\t" << (**iter_ig).getEdgeFX() << std::endl;
    }
    description << std::endl << "SUMMARY STATISTICS" << std::endl << std::endl;

    description << "Sample Size: " << m_sample_size << std::endl;

    description << "Marker Spacing: " << m_marker_spacing << std::endl;
    description << "Summary Generations:";
    std::vector<int>::iterator iterI;
    for(iterI = m_summary_generations.begin() ; iterI < m_summary_generations.end() ; iterI++ )
    {
        description << " " << (*iterI) ;
    }
    description << std::endl << "Summary Stats Collected:" << std::endl;
    description << "Allele Freq Clines: " << m_clines_on << std::endl;
    description << "Genotype Freq Clines: " << m_gclines_on << std::endl;
    description << "Marker PCA: " << m_pca_on << std::endl;
    description << "Junction Density: " << m_junction_density_on << std::endl;
    description << "Junction Frequency: " << m_jun_freq_on << std::endl;
    description << "Ancestry Blocks: " << m_blocks_on << std::endl;

    description.close();
}

// ITERATE THROUGH THE GENERATIONS
void Experiment::experimentRun()
{
    std::cout << "Experiment Running" << std::endl;
    // make sure the vector of summary generations are in order, and set up an iterator for them
    sort( m_summary_generations.begin(), m_summary_generations.end() );
    std::vector<int>::iterator iter_s = m_summary_generations.begin();
    Marker::s_generateKeys();

    for ( m_sim_number = 0 ; m_sim_number < m_number_of_simulations ; m_sim_number++ )
    {
        for ( g_generation = 0 ; g_generation < m_number_of_generations ; g_generation++ )
        {
            (*m_hap_hazard).newGeneration();

            //memMonitor();

            if ( g_generation == (*iter_s) )
            {
                std::cout << "Calculating Summaries..." << std::endl;
                Individual** sample = (*m_hap_hazard).sampleMetapop(m_sample_size);
                std::cout << "MetaPop Sampled!" << std::endl;
                summary(sample);
                std::cout << "Summary Stats Calculated!" << std::endl;
                iter_s++;
            }
        }
    }

    //memClean();
}

// CLEAN OUT THE MEMORY CONSUMED BY THE EXPERIMENT
void Experiment::memClean()
{
    std::cout << "Cleaning out memory" << std::endl;
    std::vector<Gamete*>::iterator iter_g;
    for(iter_g = m_gamete_pool.begin() ; iter_g < m_gamete_pool.end() ; iter_g++ )
    {
        //(**iter_g).recycle();
    }
    //Gamete::drainGameteRes();
    std::cout << "Gamete Pool cleared." << std::endl;
    delete m_hap_hazard;
    std::cout << "Metapopulation cleared." << std::endl;
    Individual::s_drainIndividualReservoir();
    std::cout << "Individual reservoir cleared." << std::endl;
    Chromosome::s_drainChromosomeReservoir();
    std::cout << "Chromosome Reservoir cleared." << std::endl;
    CNode::drainCNodeReservoir();
    std::cout << "CNode reservoir cleared." << std::endl;
    Junction::s_drainJunctionPool();
    std::cout << "Junction Pool cleared." << std::endl;
    Junction::s_drainJunctionReservoir();
    std::cout << "Junction Reservoir cleared." << std::endl;

    std::cout << "Memory cleared." << std::endl;

}

// PRINT OUT THE SIZE OF MEMORY RESERVE OBJECTS
void Experiment::memMonitor()
{
    std::cout << "MEMORY RESERVED AND USED" << std::endl;
    std::cout << "Jun_pool: " << Junction::s_junction_pool.size() << std::endl;
    std::cout << "Jun_res: " << Junction::s_junction_reservoir.size() << std::endl;
    std::cout << "Jun_recycled: " << Junction::s_num_recycled << std::endl;
    std::cout << "Num_Junctions: " << Junction::s_num_junctions << std::endl;
    std::cout << "CNode_res: " << CNode::s_cnode_reservoir.size() << std::endl;
    std::cout << "Num_CNode: " << CNode::s_num_cnode << std::endl;
    std::cout << "Chr_res: " << Chromosome::s_chromosome_reservoir.size() << std::endl;
    std::cout << "Ind_res: " << Individual::s_individual_reservoir.size() << std::endl;
    std::cout << "Num_Gametes: " << Gamete::s_num_gametes << std::endl;
}

// GENERATE SUMMARIES OF THE DATA
void Experiment::summary(Individual** samp)
{
    genotypeSample(samp);

    if(m_clines_on)
    {
        std::string e = "clines";

        const char* clines_file = makeFileName(e);
        std::cout << clines_file << std::endl;
        Marker::s_cline_file.open(clines_file);

        Marker::s_cline_file << "chr,pos" ;

        for( int d = 0 ; d < m_number_of_demes ; d++ )
        {
            for( int a = 0 ; a < m_number_of_ancestries ; a++ )
            {
                Marker::s_cline_file << ",D" << d << "_A" << a ;
            }
        }
        Marker::s_cline_file << std::endl;


        std::vector<Marker*>::iterator iter_m;
        for(iter_m = m_markers.begin() ; iter_m < m_markers.end() ; iter_m++ )
        {
            (**iter_m).clinesMarker();
        }
        Marker::s_cline_file.close();
    }

    if(m_gclines_on)
    {
        std::string e = "gclines";
        const char* gClinesFile = makeFileName(e);
        Marker::s_g_cline_file.open(gClinesFile);

        Marker::s_g_cline_file << "chr,pos" ;

        std::vector<std::string>::iterator iter_k;
        for( iter_k = Marker::s_genotype_keys.begin() ; iter_k < Marker::s_genotype_keys.end() ; iter_k++ )
        {
            Marker::s_g_cline_file << ",G_" << (*iter_k);
        }

        Marker::s_g_cline_file << std::endl;

        std::vector<Marker*>::iterator iter_m;
        for(iter_m = m_markers.begin() ; iter_m < m_markers.end() ; iter_m++ )
        {
            (**iter_m).gclinesMarker();
        }
        Marker::s_g_cline_file.close();
    }

    if(m_blocks_on)
    {
        blockSummary(samp);
    }


}

// GENERATE SUMMARIES OF THE DATA
void Experiment::dbSummary(Individual** samp)
{
    genotypeSample(samp);

    if(m_clines_on)
    {
        std::string e = "clines";

        const char* clines_file = makeFileName(e);
        std::cout << clines_file << std::endl;
        Marker::s_cline_file.open(clines_file);

        Marker::s_cline_file << "chr,pos" ;

        for( int d = 0 ; d < m_number_of_demes ; d++ )
        {
            for( int a = 0 ; a < m_number_of_ancestries ; a++ )
            {
                Marker::s_cline_file << ",D" << d << "_A" << a ;
            }
        }
        Marker::s_cline_file << std::endl;


        std::vector<Marker*>::iterator iter_m;
        for(iter_m = m_markers.begin() ; iter_m < m_markers.end() ; iter_m++ )
        {
            (**iter_m).clinesMarker();
        }
        Marker::s_cline_file.close();
    }

    if(m_gclines_on)
    {
        std::string e = "gclines";
        const char* gClinesFile = makeFileName(e);
        Marker::s_g_cline_file.open(gClinesFile);

        Marker::s_g_cline_file << "chr,pos" ;

        std::vector<std::string>::iterator iter_k;
        for( iter_k = Marker::s_genotype_keys.begin() ; iter_k < Marker::s_genotype_keys.end() ; iter_k++ )
        {
            Marker::s_g_cline_file << ",G_" << (*iter_k);
        }

        Marker::s_g_cline_file << std::endl;

        std::vector<Marker*>::iterator iter_m;
        for(iter_m = m_markers.begin() ; iter_m < m_markers.end() ; iter_m++ )
        {
            (**iter_m).gclinesMarker();
        }
        Marker::s_g_cline_file.close();
    }

    if(m_blocks_on)
    {
        blockSummary(samp);
    }


}

// MAKE MARKERS TO BE GENOTYPED IN THE METAPOPULATION
void Experiment::makeMarkers(double s)
{
    for(int i = 0 ; i < m_number_of_chromosomes ; i++ )
    {
        double length = 0;

        if (m_chromosome_types[i] == M || m_chromosome_types[i] == C || m_chromosome_types[i] == CP )
        {
            length = 0;
        }
        else
        {
            length = m_chromosome_lengths[i];
        }

        
        for(double j = 0 ; j <= length ; j += s)
        {
            Marker* newMarker = new Marker(i, j, m_chromosome_types[i]);

            m_markers.push_back(newMarker);
            //(*newMarker).error();
        }
    }
}

// GENOTYPE THE MARKERS USED FOR SUMMARIES OF THE EXPERIMENT
void Experiment::genotypeSample(Individual** sample)
{
    std::vector<Marker*>::iterator iter_m;
    for(iter_m = m_markers.begin() ; iter_m < m_markers.end() ; iter_m++ )
    {
        (**iter_m).setGenoToZero();
    }

    for( int d = 0 ; d < m_number_of_demes ; d++ )
    {
        for( int n = 0 ; n < m_sample_size ; n++ )
        {
            for ( iter_m = m_markers.begin() ; iter_m < m_markers.end() ; iter_m++ )
            {
                (**iter_m).g_Ind( sample[m_sample_size*d + n] );
            }
        }
    }
}

// MAKE A FILENAME FOR OUTPUT WITH EXTENSION E
const char* Experiment::makeFileName(std::string e )
{
    std::stringstream fileString;

    fileString << m_experiment_name << "_" << m_job_id << "/" << m_experiment_name << "." << m_sim_number << "." << g_generation << "." << e ;

    std::string file = fileString.str();

    return file.c_str();;
}

// MAKE A MARKER FILENAME FOR OUTPUT WITH EXTENSION E
const char* Experiment::makeMarkerFileName(std::string e )
{
    std::stringstream fileString;

    fileString << m_experiment_name << "_" << m_job_id << "/" << e << "/" << m_sim_number << "." << g_generation << "." << e ;

    std::string file = fileString.str();

    return file.c_str();
}

// MAKE A MARKER FILENAME FOR OUTPUT WITH EXTENSION E
std::string Experiment::makeBlockFileName(std::string e, int c)
{
    std::stringstream fileString; 
    
    fileString << m_experiment_name << "_" << m_job_id << "/BL_CHR" << c << "/" << m_sim_number << "." << g_generation << "." << e ;
    std::cout << "EX_589:" << fileString.str() << std::endl;
    std::string file = fileString.str();

    std::ofstream bf_list;
    bf_list.open("block_file_list.txt", std::ios_base::app);
    bf_list << file.c_str() << "," << m_marker_spacing << "," << m_chromosome_lengths[c] << std::endl;
    bf_list.close();

    return file;
}


// GET THE ANCESTRY BLOCKS FROM A SAMPLE OF THE METAPOPULATION -- print them straight to file
void Experiment::blockSummary(Individual** sample)
{
    std::string e = "bla";
    std::cout << "Block Summary: " << std::endl;
    for(int c = 0 ; c < m_number_of_chromosomes ; c++ )
    {
        std::string block_file = makeBlockFileName(e, c);

        std::cout << block_file << std::endl;
        m_os_blocks.open(block_file.c_str());
        m_os_blocks << "Deme,Individual,Haplotype,Chromosome,Type,Ancestry,Start,Start_Jun,Start_Born,End,Length" << std::endl;
        if( m_chromosome_types[c] != M && m_chromosome_types[c] != C && m_chromosome_types[c] != CP )
        {
            for( int i = 0 ; i < m_sample_size * m_number_of_demes ; i++ )
            {
                std::vector<Block> ind_blocks = (*sample[i]).getBlocks(c);

                std::vector<Block>::iterator iter_b;
                for( iter_b = ind_blocks.begin() ; iter_b < ind_blocks.end() ; iter_b++ )
                {

                    // deme, individual, chromosome, haplotype, ctype, ancestry, start, end, length
                    m_os_blocks << (*sample[i]).getLocation() ;
                    m_os_blocks << "," << (sample[i]) ;
                    m_os_blocks << "," << (*iter_b).haplotype ;
                    m_os_blocks << "," << (*iter_b).chromosome ;
                    m_os_blocks << "," << (*iter_b).cType ;
                    m_os_blocks << "," << (*iter_b).ancestry ;
                    m_os_blocks << "," << (*iter_b).start ;
                    m_os_blocks << "," << (*iter_b).startJunction ;
                    m_os_blocks << "," << (*iter_b).startGen ;
                    m_os_blocks << "," << (*iter_b).end ;
                    m_os_blocks << "," << (*iter_b).end - (*iter_b).start << std::endl;

                }
            }
        }
        m_os_blocks.close();
    }
    std::cout << "Block Summary Done!" << std::endl;
}

// SETUP AND INTIALIZE THE OBJECTS NEEDED FOR THE EXPERIMENT BASED ON THE INPUT PARAMETERS
void Experiment::experimentSetup()
{
    // Make the gamete pool for the initial generations/new individuals etc.
    std::cout << "Checking mating system and building gamete pool... " << std::endl;
    switch(m_mating_system)
    {
        case(HH):
            m_gamete_pool = HHgametePool(m_number_of_ancestries, m_chromosome_lengths, m_chromosome_types); //(int anc, vector<double> lengths, vector<ChrType> types, double p = 0.001)
            break;
        case(XY):
            m_gamete_pool = XYgametePool(m_number_of_ancestries, m_chromosome_lengths, m_chromosome_types, m_par_length);
            break;
        case(ZW):
            m_gamete_pool = ZWgametePool(m_number_of_ancestries, m_chromosome_lengths, m_chromosome_types, m_par_length);
            break;
        default:
            std::cout << "The gamete pool could not be made. The mating system was not properly specified." << std::endl;
            exit(2);
    }

    // set gamete pool individual static variable
    Individual::s_ind_gamete_pool = m_gamete_pool;

    std::vector<Gamete*>::iterator iter_gp;
    for(iter_gp = m_gamete_pool.begin() ; iter_gp < m_gamete_pool.end() ; iter_gp++ )
    {
        (**iter_gp).displayGamete();
    }

    std::cout << "... Done!" << std::endl;

    std::cout << std::endl;
    std::cout << "Building genotype-phenotype landscapes... ";
    //Landscape( double a, double exp_value, int na, vector<Gene*> l, vector<IntGraph*> iG );
    // Build a landscape for each phenotype
    for ( int i = 0 ; i < m_number_of_phenotypes ; i++ )
    {
        std::vector<Gene*>::iterator iter_g;
        std::vector<Gene*> landscapeGenes;
        for( iter_g = m_genes.begin() ; iter_g < m_genes.end() ; iter_g++ )
        {
            int j = (**iter_g).getPhenotype();
            if( j == i )
            {
                landscapeGenes.push_back(*iter_g);
            }
        }

        std::vector<IntGraph*>::iterator iter_i;
        std::vector<IntGraph*> landscapeNetwork;
        for( iter_i = m_network.begin() ; iter_i < m_network.end() ; iter_i++ )
        {
            int j = (**iter_i).getPhenotype();
            if( j == i )
            {
                landscapeNetwork.push_back(*iter_i);
            }
        }

        Landscape* newLandscape = new Landscape(i, m_exp_phen_values[i], m_number_of_ancestries, landscapeGenes, landscapeNetwork);
        (*newLandscape).printLandscape();
        m_phen_landscapes.push_back(newLandscape);
    }

    Individual::s_landscapes = m_phen_landscapes;

    std::cout << "Done!" << std::endl;

    // BUILD THE METAPOP
    // MetaPop(int d, vector<int> s, vector<double> a, vector<Gamete*> g, MatingSys ms, vector<double> mr)

    std::cout << "Setting up the Metapopulation...";
    m_hap_hazard = new MetaPop( m_number_of_demes, m_deme_sizes, m_ancestry_frequencies, m_gamete_pool, m_mating_system, m_migration_rates);
    std::cout << " Done!" << std::endl;

    std::cout << "Setting up output parameters...";
    //Make markers
    makeMarkers(m_marker_spacing);
    //Make the junction density table
    double genomeLength = 0;
    std::vector<double>::iterator iterL;
    for( iterL = m_chromosome_lengths.begin() ; iterL < m_chromosome_lengths.end() ; iterL++ )
    {
        genomeLength += (*iterL);
    }
    std::cout << " Done!" << std::endl;
}

// ENTER PARAMETERS FROM A FILE
void Experiment::fileInput()
{
    std::ifstream input;
    input.open(m_INPUTFILE);
    char inpBuffer[30];
    int buf = 30;

    //std::cout << "What would you like to name your experiment? <text> (Your results will be stored in a directory of the same name.): ";
    input.getline(inpBuffer, buf);
    m_experiment_name = inpBuffer;

    //std::cout << "How many simulations will be run? <integer>: ";
    input.getline(inpBuffer, buf);
    m_number_of_simulations = atoi(inpBuffer);

    //std::cout << "How many generations should be in each simulation? <integer>: ";
    input.getline(inpBuffer, buf);
    m_number_of_generations = atoi(inpBuffer);

    //std::cout << "How many demes (sub-poplations) will be used in your simulations? <integer>: ";
    input.getline(inpBuffer, buf);
    m_number_of_demes = atoi(inpBuffer);

    //std::cout << "Now enter the list of populations sizes for each deme in order. <integer>" << std::endl;
    for (int i = 0 ; i < m_number_of_demes ; i++ )
    {
        //std::cout << "Size of Deme " << i << ": ";
        input.getline(inpBuffer, buf);
        int n = atoi(inpBuffer);
        m_deme_sizes.push_back(n);
    }

    //std::cout << "How many ancestries are in the metapopulation? <integer>:";
    input.getline(inpBuffer, buf);
    m_number_of_ancestries = atoi(inpBuffer);

    //std::cout << "Please enter the frequency of each ancestry, 0 through " << (m_number_of_ancestries - 1) << ", enter its frequency in each deme. <decimal number, double>" << std::endl;
    //std::cout << "Ancestry frequencies in each deme MUST sum to 1.0" << std::endl;

    for( int i = 0 ; i < m_number_of_demes ; i++ )
    {
        for( int j = 0 ; j < m_number_of_ancestries ; j++ )
        {
            //std::cout << "Frequency of Ancestry " << j << "in Deme " << i << ": ";
            input.getline(inpBuffer, buf);
            double a = strtod(inpBuffer, NULL);
            m_ancestry_frequencies.push_back(a);
            //std::cout << i << " " << j << " " << a << std::endl;
        }
    }

    //std::cout << "What type of mating/sex deterination system does the population have? <integer>" << std::endl;
    //std::cout << "<0> = Hermaphrodites, <1> = XX females and XY males, <2> = WZ females and ZZ males" << std::endl;
    input.getline(inpBuffer, buf);
    int MS = atoi(inpBuffer);
    m_mating_system = MatingSys(MS);

    //std::cout << "Now enter the migrations rates between demes. <decimal/double>" << std::endl;
    //std::cout << "Start and end with the migration rate with source populations." << std::endl;

    for( int i = 0 ; i <= m_number_of_demes ; i++ )
    {
        if (i == 0)
        {
            //std::cout << "Source A and Deme 0: ";
            input.getline(inpBuffer, buf);
            double m = strtod(inpBuffer, NULL);
            m_migration_rates.push_back(m);
        }
        else if ( i == m_number_of_demes )
        {
            //std::cout << "Source B and Deme " << i << ": ";
            input.getline(inpBuffer, buf);
            double m = strtod(inpBuffer, NULL);
            m_migration_rates.push_back(m);
        }
        else
        {
            //std::cout << "Deme " << i-1 << " and Deme " << i-1 << ": ";
            input.getline(inpBuffer, buf);
            double m = strtod(inpBuffer, NULL);
            m_migration_rates.push_back(m);
        }
    }

    //std::cout << "Next enter the genomic parameters." << std::endl;
    //std::cout << "How many chromosomes will the genome contain (in a haplotype)? <integer> (including up to one sex chromosome, and one cytoplasmic factor): ";
    input.getline(inpBuffer, buf);
    m_number_of_chromosomes = atoi(inpBuffer);

    //std::cout << "Which model of recombination would you like to use?" << std::endl;
    //std::cout << "Options: <0> = one crossover per chromosome, <1> = Poisson recombination without interference, <2> = Gamma distributed interference: ";
    input.getline(inpBuffer, buf);
    m_interference_option = atoi(inpBuffer);

    //std::cout << "What types of chromosomes will be used?" << std::endl;
    //std::cout << "Options (case sensitive):<0> = autosomal pair, <2> = XY sex chromosome pair, <4> = ZW sex chromosome pair,  <5> = mitochondrion, <6> = chloroplast, <7> = paternally inherited cytoplasm " << std::endl;
    //std::cout << "IMPORTANT NOTE: the sexchromosome pair, if used MUST be listed first, followed by any number of autosomes, and the cytoplams MUST be listed last. " << std::endl;
    //std::cout << "Not following this format may cause the program to crash." << std::endl;

    //std::cout << "Number of Chromosomes is " << m_number_of_chromosomes << std::endl;

    for( int i = 0 ; i < m_number_of_chromosomes ; i++ )
    {
       //std::cout << "Enter the type for chromosome " << i << ": ";
        input.getline(inpBuffer, buf);
        int ct = atoi(inpBuffer);
        ChrType t = ChrType(ct);
        m_chromosome_types.push_back(t);
    }

    //std::cout << "What are the genetic lengths (in Morgans) of these chromosomes? <decimal/double>" << std::endl;
    //std::cout << "Follow the same order as the types. Ignore the cytoplasm, it will automatically be set to zero." << std::endl;

    for( int i = 0 ; i < m_number_of_chromosomes ; i++ )
    {
        //std::cout << "Enter the genetic length of chromosome " << i << ": " << std::endl;
        input.getline(inpBuffer, buf);
        double l = strtod(inpBuffer, NULL);
        m_chromosome_lengths.push_back(l);
    }

    //std::cout << "What is the genetic lenght of the pseudoautosomal region? <decimal/double>" << std::endl;
    //std::cout << "If you are not using sex chromosomes this paramter will be ignored. If there are sex chromosomes, but you wish to ignore it enter <0>." << std::endl;
    input.getline(inpBuffer, buf);
    m_par_length = strtod(inpBuffer, NULL);

    //std::cout << "Now the genetic parameters and phenotypes." << std::endl;

    //std::cout << "How many phenotypes? <integer>" << std::endl;
    //std::cout << "Phenotypes 1, 2, and 3 will automatically be reproductive, developmental, and environmental fitness respectively." << std::endl;
    //std::cout << "Please enter a 3 for the defaults, or 4 or more if you wish to track other phenotypes. You can ignore the fitnesses in your gene-phenotype assignments if you wish." << std::endl;
    input.getline(inpBuffer, buf);
    m_number_of_phenotypes = atoi(inpBuffer);

    for (int i = 0 ; i < 3 ; i++ )
    {
        m_exp_phen_values.push_back(1);
    }

    if( m_number_of_phenotypes > 3 )
    {
        //std::cout << "The expected values of the fitnesses will be automatically set to 1.0." << std::endl;
        //std::cout << "Please enter the expected values of the extra phenotypes you have specified. <decimal/double>" << std::endl;
        for( int i = 3 ; i < m_number_of_phenotypes ; i++ )
        {
            //std::cout << "The expected value for phenotype " << i << " = ";
            input.getline(inpBuffer, buf);
            double epv = strtod(inpBuffer, NULL);
            m_exp_phen_values.push_back(epv);
        }
    }

    // INITIALIZE STATIC CLASS VARIABLES
    // Random Number Generator Variables
    MetaPop::s_rand_num = s_rand_num;
    Deme::s_rand_num = s_rand_num;
    Individual::s_rand_num = s_rand_num;
    Chromosome::s_rand_num = s_rand_num;

    // Set the par boundary and junction density window size
    Chromosome::s_par_boundary = m_par_length;
    Chromosome::s_num_anc = m_number_of_ancestries;

    // Set the number of phenotypes individuals can hold
    Individual::s_num_of_phenotypes = m_number_of_phenotypes;
    Individual::s_num_chr = m_number_of_chromosomes;
    Individual::s_interference_opt = m_interference_option;
    Individual::s_gp_type = m_mating_system;
    Individual::s_min_block_size = 0.00000001;

    // Set the number of gametes to zero
    Gamete::s_num_gametes = 0;

    // Set deme static variables
    Deme::s_mat_sys = m_mating_system;
    Deme::s_num_anc = m_number_of_ancestries;

    //Set gene static variables
    Gene::s_num_anc = m_number_of_ancestries;
    Gene::s_num_chr = m_number_of_chromosomes;

    // Set intgraph static variables
    IntGraph::s_num_ancestries = m_number_of_ancestries;

    // Set the number of demes for the junction frequency array
    Junction::s_num_locations = m_number_of_demes;
    Junction::s_num_junctions = 0;

    // CNode static variables
    CNode::s_num_cnode = 0;

    // Set Marker static functions
    Marker::s_num_demes = m_number_of_demes;
    Marker::s_num_ancestries = m_number_of_ancestries;
    Marker::s_mark_sample_size = m_sample_size;
    Marker::s_generateKeys();

    //std::cout << "How many genes do you need in your simulations?" << std::endl;
    input.getline(inpBuffer, buf);
    int numGenes = atoi(inpBuffer);

    for( int g = 0 ; g < numGenes ; g++ )
    {
        //std::cout << "Which chromosome is gene " << g << " on? <integer>: ";
        input.getline(inpBuffer, buf);
        int chr = atoi(inpBuffer);

        //std::cout << "And which position? <decimal/double>: ";
        input.getline(inpBuffer, buf);
        double pos = strtod(inpBuffer, NULL);

        //std::cout << "Which Phenotype does it affect? <0> = reproductive fitness, <1> = developmental fitness, <2> = environmental fitness, <3> = other phenotype: ";
        input.getline(inpBuffer, buf);
        int pheno = atoi(inpBuffer);

        //std::cout << "And what is the additive effect of each ancestry on the phenotype? <decimal/double> (must be 0 <= a <= 1 for fitnesses): ";
        std::vector<double> addFX;

        int anc_multiply = 1;
        if( chr == 0)
        {
            anc_multiply = 2;
        }

        for (int j = 0 ; j < m_number_of_ancestries * anc_multiply ; j++ )
        {
            input.getline(inpBuffer, buf);
            double aE = strtod(inpBuffer, NULL);
            addFX.push_back(aE);
        }

        Gene* newGene = new Gene(chr, pos, addFX, pheno);
        m_genes.push_back(newGene);
    }

    //std::cout << "Are there any interactions between pairs of genes? If so how many? <integer>: ";
    input.getline(inpBuffer, buf);
    int numInt = atoi(inpBuffer);

    for( int i = 0 ; i < numInt ; i++ )
    {
        //std::cout << "Which two genes participate in interaction " << i << "?" << std::endl;
        //std::cout << "Specify the indices (0-n) of the genes in the list you just made. <integer>: ";
        //std::cout << "The first gene is: ";
        input.getline(inpBuffer, buf);
        int iOne = atoi(inpBuffer);
        //std::cout << "The second gene is: ";
        input.getline(inpBuffer, buf);
        int iTwo = atoi(inpBuffer);
        //std::cout << "Next, specify which ancestries of the genes interact. <integer>" << std::endl;
        //std::cout << "The interacting ancestry at the first gene is: " << std::endl;
        input.getline(inpBuffer, buf);
        int aOne = atoi(inpBuffer);
        //std::cout << "The interacting ancestry at the second gene is: " << std::endl;
        input.getline(inpBuffer, buf);
        int aTwo = atoi(inpBuffer);
        //std::cout << "What is the maximum phenotypic effect of this interaction? <decimal/double>: " << std::endl;
        input.getline(inpBuffer, buf);
        double pE = strtod(inpBuffer, NULL);
        //std::cout << "Which phenotype does the interaction affect? <integer> " << std::endl;
        input.getline(inpBuffer, buf);
        int iP = atoi(inpBuffer);
        //std::cout << "What model of epistasis does it follow? <integer>, <0> = recessive, <1> = dominant, <2> = dom-rec, <3> = additive: " << std::endl;
        input.getline(inpBuffer, buf);
        int mod = atoi(inpBuffer);

        IntGraph* newInteraction = new IntGraph( iP, m_genes[iOne], m_genes[iTwo], aOne, aTwo, pE, mod);
        m_network.push_back(newInteraction);
    }

    // Output Files and Summary Statistics
    //std::cout << "Now we'll decide how to get the data you need." << std::endl;
    //std::cout << "How many samples should be collected from each deme? <integer>" << std::endl;
    input.getline(inpBuffer, buf);
    m_sample_size = atoi(inpBuffer);
    //std::cout << "How big should the windows for junction densities be? <decimal/double>" << std::endl;
    //input.getline(inpBuffer, buf);
    //windowSize = strtod(inpBuffer, NULL);

    //std::cout << "How far apart should markers be spaced? <double>" << std::endl;
    input.getline(inpBuffer, buf);
    m_marker_spacing = strtod(inpBuffer, NULL);

    //std::cout << "Please enter the integers identifying the generations at which you would like to collect summaries. <integer>:" << std::endl;
    //std::cout << "Enter a negative number or a number beyond the range of the total number of generations to quit." << std::endl;

    int summaryInt = 0;

    while( summaryInt >= 0 && summaryInt <= m_number_of_generations )
    {
        input.getline(inpBuffer, buf);
        summaryInt = atoi(inpBuffer);
        m_summary_generations.push_back(summaryInt);
        //std::cout << "Ex_FI()_990 " << summaryInt << " " << m_number_of_generations << std::endl;
    }

    //std::cout << "Finally, which data would you like to collect? Enter 0 to turn data off, enter 1 to turn it on." << std::endl;
    //std::cout << "Geographic allele frequency clines?: ";
    input.getline(inpBuffer, buf);
    m_clines_on = atoi(inpBuffer);

    //std::cout << "Geographic genotype frequency clines?: ";
    input.getline(inpBuffer, buf);
    m_gclines_on = atoi(inpBuffer);

    //std::cout << "Marker PCA?: " ;
    input.getline(inpBuffer, buf);
    m_pca_on = atoi(inpBuffer);

    //std::cout << "Junction Densities?: ";
    input.getline(inpBuffer, buf);
    m_junction_density_on = atoi(inpBuffer);

    //std::cout << "Junction Frequencies?: ";
    input.getline(inpBuffer, buf);
    m_jun_freq_on = atoi(inpBuffer);

    //std::cout << "Ancestry Blocks?: ";
    input.getline(inpBuffer, buf);
    m_blocks_on = atoi(inpBuffer);
}

#endif // EXPERIMENT_H
