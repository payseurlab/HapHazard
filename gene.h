#ifndef GENE_H
#define GENE_H

#include <vector>
#include <algorithm>
#include <cmath>

#include "chromosome.h"
#include "jungen_acc.cpp"

/*
    GENE

    Genes store information that defines a gene, and can genotype an individual.
*/

class Gene
{
    public:
        Gene(int chr, double pos, std::vector<double> ae, int phen);


        // Getters
        int getChr() { return m_chromosome; }
        double getPos() { return m_position; }
        double getAddEffect(int i) { return m_add_effects[i]; }
        std::vector<double> getAddFX() { return m_add_effects; }
        int getDominance(int i) { return m_dominance[i]; }
        int getPhenotype() { return m_phenotype; }
        std::vector<std::vector<int> > getGenotypes() { return m_genotypes; }
        std::vector<int> getGenotype(int i) { return m_genotypes[i]; }

        // Setters


        // Biological Functions
        std::vector<int> genotype( std::vector<Chromosome*> );

        // STATIC VARIABLES
        static int s_num_anc;
        static int s_num_chr;

    protected:
        int m_chromosome;                   // the chromosome the gene is on
        double m_position;                  // the positon of the chromosome
        std::vector<double> m_add_effects;       // the additive effects of the gene on its phenotype
        std::vector<int> m_dominance;            // the dominance of a gene's effects
        int m_phenotype;                    // the phenotype the gene affects
        std::vector<std::vector<int> > m_genotypes;   // the list of possible genotypes the gene can have

};

// GENE CONSTRUCTOR
Gene::Gene(int chr, double pos, std::vector<double> ae, int p = 0):
    m_chromosome(chr), m_position(pos), m_add_effects(ae), m_phenotype(p)
{
    int num_of_genotypes = 0;
    int anc_multi = 1;
    if( m_chromosome == 0 )
    {
        num_of_genotypes = pow( 2 * s_num_anc, 2);
        anc_multi = 2;
    }
    else if( m_chromosome == s_num_chr - 1)
    {
        num_of_genotypes = s_num_anc;
    }
    else
    {
        num_of_genotypes = pow( s_num_anc, 2);
    }

    for( int i = 0 ; i < num_of_genotypes ; i++ )
    {
        std::vector<int> new_genotype = change_base(i, s_num_anc*anc_multi, 2);
        std::vector<int>::iterator iter_i;
        for(iter_i = new_genotype.begin() ; iter_i < new_genotype.end() ; iter_i++ )
        {

        }
        m_genotypes.push_back(new_genotype);
    }

}


// GET THE GENOTYPE FOR A GENE IN AN INDIVIDUAL
std::vector<int> Gene::genotype(std::vector<Chromosome*> g)
{
    std::vector<int> genotype;
    int x = m_chromosome * 2;

    genotype.push_back( (*g[x]).positionAnc(m_position) );
    genotype.push_back( (*g[x+1]).positionAnc(m_position) );

    sort( genotype.begin(), genotype.end() );

    return genotype;
}


#endif // GENE_H
