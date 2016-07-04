#ifndef MARKER_H
#define MARKER_H

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>

#include "chromosome.h"
#include "individual.h"
#include "deme.h"
#include "metaPop.h"

/*
    MARKER


*/

extern std::ofstream g_check_file;

class Marker
{
    public:
        // CONSTRUCTOR PROTOTYPE
        Marker(int c, double p, ChrType t);

        // Getters
        double getPosition() { return m_mark_position; }
        int getChromosome() { return m_mark_chromosome; }
        //vector<double> getDemeFreq() { return demeFreqs; }
        double getMPFreq() {return m_mp_freq; }
        int* getLocationCounts() { return m_location_counts; }
        std::map<std::string,int>* getGenotypeCounts() { return &m_genotype_counts; }

        // Setters
        void setPosition ( double p ) { m_mark_position = p; }
        void setChromosome ( int c ) { m_mark_chromosome = c; }
        //void setDemeFreqs ( vector<double> f ) { demeFreqs = f; }

        // Genotyping functions
        int g_Chr(Chromosome* c);                    // genotype the position on a single chromosome
        void g_Ind(Individual* i);           // genotype the position on both chromosomes of an individual
        double g_Deme(Deme* d);                     // determine the ancestry frequency in given deme
        void g_Demes(MetaPop* mp);                   // determine the ancestry frequencies in all demes of a metapopulation -- returns it to demeFreq
        void g_Meta(MetaPop* mp);                  // determine the ancestry frequency in a metapopulation
        void gVar_Meta();               // determine the variance in ancestry freq among demes in a metapopulation
        void setGenoToZero();

        // static members/functions
        static int s_num_ancestries;
        static int s_num_demes;
        static int s_mark_sample_size;
        static std::vector<std::string> s_generateKeys();
        static std::vector<std::string> s_genotype_keys;
        static std::ofstream s_cline_file;
        static std::ofstream s_g_cline_file;

        // Data Management functions
        void clinesMarker();
        void gclinesMarker();
        void error();

    private:
        int m_mark_chromosome;            // the chromosome the marker is on
        double m_mark_position;           // the physical/genomic position of the marker
        double m_mp_freq;                  // the marker ancestry frequency in the meta-population
        double m_mp_var_freq;               // the variance in ancestry frequency among demes
        int* m_location_counts;
        std::map<std::string,int> m_genotype_counts;
};

// MARKER CONSTRUCTOR
Marker::Marker(int c = -1, double p = -1, ChrType t = A):
    m_mark_chromosome(c), m_mark_position(p)
{
    if( t == A || t == C || t == M || t == CP )
    {
        m_location_counts = new int[s_num_demes * s_num_ancestries];
        for( int i = 0 ; i < s_num_demes * s_num_ancestries ; i++ )
        {
            m_location_counts[i] = 0;
        }
    }
    else
    {
        m_location_counts = new int[s_num_demes * s_num_ancestries * 2 ];
        for( int i = 0 ; i < s_num_demes * s_num_ancestries * 2 ; i++ )
        {
            m_location_counts[i] = 0;
        }
    }

    std::map<std::string,int> m_genotype_counts;

    std::vector<std::string>::iterator iter_k;
    for( iter_k = s_genotype_keys.begin() ; iter_k < s_genotype_keys.end() ; iter_k++ )
    {
        m_genotype_counts[*iter_k] = 0;
    }
}

// SET ALL THE GENOTYPE COUNTS TO ZERO
void Marker::setGenoToZero()
{
    std::map<std::string,int>::iterator iter_g;
    for( iter_g = m_genotype_counts.begin() ; iter_g != m_genotype_counts.end() ; iter_g++ )
    {
        (*iter_g).second = 0;
    }

    for( int i = 0 ; i < s_num_demes * s_num_ancestries ; i++ )
    {
        m_location_counts[i] = 0;
    }
}

// PRINT THE MARKER TO ERROR CHECK
void Marker::error()
{
    g_check_file << "Marker: " << m_mark_chromosome << " " << m_mark_position << std::endl;
}

// PRINT THE MARKERS TO A CLINES OUTPUT FILE
void Marker::gclinesMarker()
{
    Marker::s_g_cline_file << m_mark_chromosome << "," << m_mark_position;

    std::vector<std::string>::iterator iter_k;
    for( iter_k = Marker::s_genotype_keys.begin() ; iter_k < Marker::s_genotype_keys.end() ; iter_k++ )
    {
        s_g_cline_file << "," << m_genotype_counts[*iter_k];
    }

    s_g_cline_file << std::endl;
}

// PRINT THE MARKERS TO A GENEOTYPE CLINES OUTPUT FILE
void Marker::clinesMarker()
{
    Marker::s_cline_file << m_mark_chromosome << "," << m_mark_position;

    for( int i = 0 ; i < s_num_demes*s_num_ancestries; i++ )
    {

        Marker::s_cline_file << "," << m_location_counts[i];
    }

    Marker::s_cline_file << std::endl;

}
// GENOTYPING FUNCTIONS

// CREATE A LIST OF MARKER GENOTYPE KEYS
std::vector<std::string> Marker::s_generateKeys()
{
    for( int t = 0 ; t <= s_num_ancestries ; t++ )
    {
        for( int d = 0 ; d < s_num_demes ; d++ )
        {
            for( int g = 0 ; g < s_num_ancestries ; g++ )
            {
                for( int h = g ; h < s_num_ancestries ; h++ )
                {
                    std::stringstream key;
                    if( t < s_num_ancestries )
                    {
                        key << d << g << h << t;
                    }
                    else
                    {
                        key << d << g << h;
                    }
                    s_genotype_keys.push_back( key.str() );
                }
            }
        }
    }
    return s_genotype_keys;
}

// GENOTYPE THE MARKER ON A SINGLE CHROMOSOME
int Marker::g_Chr(Chromosome* c)
{
    int anc = (*c).positionAnc(m_mark_position);
    return anc;
}

// GENOTYPE THE POSITION ON BOTH CHROMOSOMES OF AN INDIVIDUAL
void Marker::g_Ind(Individual* i)
{
    Chromosome* chromosome_one = (*i).getChr(m_mark_chromosome * 2);
    Chromosome* chromosome_two = 0;
    int location = (*i).getLocation() ;

    // check the types of chromosomes to be genotyped
    // if the chromosome is a cytoplasmic genome, set the chromosome_two pointer to NULL
    if( (*chromosome_one).getType() == A || (*chromosome_one).getType() == X || (*chromosome_one).getType() == Z )
    {
        chromosome_two = (*i).getChr(m_mark_chromosome * 2 + 1);
    }


    // Genotype the marker on the chromosome(s)
    // if it is in the cytoplasm, genotype one chromosome and score it as homozygous
    if( (*chromosome_one).getType() == M || (*chromosome_one).getType() == C || (*chromosome_one).getType() == CP )
    {
        int geno_one = (*chromosome_one).positionAnc(m_mark_position);
        m_location_counts[ location * s_num_ancestries + geno_one ] += 1;
        std::stringstream genotype;
        genotype << location << geno_one << geno_one << 0;
        m_genotype_counts[genotype.str()] += 1;
    }
    // otherwise, genotype it as normally, the chromosome will take care of genotyping beyond the pseudoautosomal region
    else
    {
        int geno_one = (*chromosome_one).positionAnc(m_mark_position);
        int geno_two = (*chromosome_two).positionAnc(m_mark_position);

        if( (*chromosome_one).getType() != Y && (*chromosome_one).getType() != W )
        {
            m_location_counts[ location * s_num_ancestries + geno_one ] += 1;
        }

        if( (*chromosome_two).getType() != Y && (*chromosome_two).getType() != W )
        {
            m_location_counts[ location * s_num_ancestries + geno_two ] += 1;
        }

        std::stringstream genotype;

        if( (*chromosome_two).getType() == Y || (*chromosome_two).getType() == W )
        {
            if( geno_one < geno_two )
            {
                genotype << location << geno_one << geno_two << geno_one;
            }
            else
            {
                genotype << location << geno_two << geno_one << geno_one;
            }
        }
        else
        {
            if( geno_one < geno_two )
            {
                genotype << location << geno_one << geno_two;
            }
            else
            {
                genotype << location << geno_two << geno_one;
            }
        }
        m_genotype_counts[ genotype.str() ] += 1;

    }
}

#endif // MARKER_H
