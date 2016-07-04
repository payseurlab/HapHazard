#ifndef METAPOP_H
#define METAPOP_H

#include <vector>
#include <algorithm>

#include "deme.h"

/*
    METAPOPULATION

    The metapopluation is a collection of demes amond which migration occurs.
    The metapopulation class defines the relationships among demes and regulates
    the migration among them.

    Currently, the metpopulation models a stepping stone model
*/

class MetaPop
{
    public:
        MetaPop(int d, std::vector<int> s, std::vector<double> a, std::vector<Gamete*> g, MatingSys ms, std::vector<double> mr );
        ~MetaPop();

        // setters
        void setNumDemes(int d) { m_num_demes = d; }

        // getters
        int getNumDemes() { return m_num_demes; }
        std::vector<Deme*> getDemes() { return m_demes; }
        Deme* getDeme(int i) { return m_demes[i]; }

        // Hybrid zone functions
        void stepStone();
        void migrate( Deme* d, Deme* e, double m);
        void migrateSource( Deme* d, double m , int s);
        void newGeneration();
        void longRangeDisp(Deme* d);


        // Data and summaries
        void junByWindow(double w, const char* f);
        Individual** sampleMetapop(int n);

        // Static variables
        static gsl_rng* s_rand_num;


    private:
        int m_num_demes;                    // the number of demes in the metapopulation
        std::vector<Deme*> m_demes;              // the collection of pointers to the demes
        std::vector<int> m_pop_sizes;            // the sizes of the demes
        std::vector<double> m_anc_proportions;   // the ancestry proportions of the demes
        std::vector<Gamete*> m_gamete_pool;      //
        std::vector<double> m_migration_rates;
        std::vector<double> m_source_a;
        std::vector<double> m_source_b;

};

// META POPULATION CONSTRUCTOR
MetaPop::MetaPop(int d, std::vector<int> s, std::vector<double> a, std::vector<Gamete*> g, MatingSys ms, std::vector<double> mr):
    m_num_demes(d), m_pop_sizes(s), m_anc_proportions(a), m_gamete_pool(g), m_migration_rates(mr)
{
    int x = m_anc_proportions.size() / m_num_demes; // this tells how much to increment by over the ancestry proportion vector
    int j = 0;  // this will be passed to each deme to assign it a number, then incremented
    std::vector<double>::iterator iter_a;

    // iterate over the ancestry proportions
    for( iter_a = m_anc_proportions.begin() ; iter_a < m_anc_proportions.end() - x + 1 ; iter_a+=x )
    {
        std::vector<double> deme_ancs;
        deme_ancs.assign(iter_a, iter_a + x); // the ancestry proportions or the currnet deme in the iteration

        Deme* d = new Deme(m_pop_sizes[j], deme_ancs, j, m_gamete_pool);

        m_demes.push_back(d);
        j++;
    }
}

// META POPULATION DESTRUCTOR
MetaPop::~MetaPop()
{
    // when the metapopulation is destroyed, erase each deme
    std::vector<Deme*>::iterator iter_d;
    for ( iter_d = m_demes.begin() ; iter_d < m_demes.end() ; iter_d++ )
    {
        delete *iter_d;
    }
}

// MIGRATE INDIVIDUALS AMONG DEMES IN THE WHOLE HYBRID ZONE AND SOURCE POPULATIONS -- STEPPING STONE MODEL
void MetaPop::stepStone()
{
    std::vector<Deme*>::iterator iter_d;
    std::vector<Deme*>::iterator iter_e;

    std::vector<double>::iterator iter_mig = m_migration_rates.begin();

    // iterate over the demes
    for ( iter_d = m_demes.begin() ; iter_d <= m_demes.end() ; iter_d++)
    {

        if ( iter_d == m_demes.begin() )    // when at the beginning, migrate from the source
        {
            migrateSource( (*iter_d), (*iter_mig), 0);

            iter_mig++;
        }
        else if ( iter_d == m_demes.end() ) // at the end, migrat from the other source
        {
            iter_e = iter_d - 1;
            migrateSource( (*iter_e), (*iter_mig), 1);
        }
        else // in the middle, migrate between two adjacent demes
        {

            iter_e = iter_d-1;
            migrate( *iter_e, *iter_d, (*iter_mig));
            iter_mig++;
        }
    }

}

// EXCHANGE RANDOM INDIVIUDALS BETWEEN TWO DEMES -- MIGRATION
void MetaPop::migrate( Deme* d, Deme* e, double m)
{
    double mu = (*d).getSize() * m / 2; // the number of migrants
    int migrants;

    if ( mu <= 1 ) // if the migration rate times the population size < 1, choose randomly from a poisson
    {
        mu = gsl_ran_poisson(s_rand_num, mu);
    }

    migrants = mu;

    // choose mu random immigrant and emmigrants pairs
    for ( int i  = 0 ; i < migrants ; i++ )
    {
        int x;
        int y;

        do
        {
            x = gsl_rng_uniform(s_rand_num) * (*d).getSize();
            //std::cout << "MP_150: " << x << "\t" << (*((*d).getMembers()[x])).p_migrated << std::endl;
        } while ( (*((*d).getMembers()[x])).p_migrated ); // make sure the member chosen hasn't already migrated

        do
        {
            y = gsl_rng_uniform(s_rand_num) * (*e).getSize();
            //std::cout << "MP_156: " << y << "\t" << (*((*d).getMembers()[y])).p_migrated << std::endl;
        } while ( (*((*d).getMembers()[y])).p_migrated ); // make sure the member chosen hasn't already migrated

        // swap the immigrant and emigrant
        Individual* temp_ind = (*d).getMembers()[x];
        (*d).immigrate(x, (*e).getMembers()[y]);
        (*e).immigrate(y, temp_ind);


    }

}

// IMMIGRATE AN INDIVIDUAL FROM A SOURCE POPULATION
void MetaPop::migrateSource( Deme* d, double m, int s)
{
    int migrants = (*d).getSize() * m / 2;

    for ( int i  = 0 ; i < migrants ; i++ )
    {
        int x;

        // get an emmigrant from the deme that hasn't already migrated
        do
        {
            x = gsl_rng_uniform(s_rand_num) * (*d).getSize();


        } while ( ( *( (*d).getMembers()[x] ) ).p_migrated );

        // randomly choose which sex the new individuals will be
        int sex_chr = gsl_ran_bernoulli(s_rand_num, 0.5);
        int G1 = 2 * s + 1;
        int G2 = 2 * s + 1 - sex_chr;

        // get gametes for the new individual
        Gamete* mig_chr_one = (*(Individual::s_ind_gamete_pool[G1])).copyGamete( (*d).getLocation() );
        Gamete* mig_chr_two = (*(Individual::s_ind_gamete_pool[G2])).copyGamete( (*d).getLocation() );
        Individual* migrant = new ( Individual::newAddyAssign() ) Individual(mig_chr_one, mig_chr_two, (*( (*d).getLocation() )), true );

        // recycle the gametes when done copying them
        (*mig_chr_one).recycle();
        (*mig_chr_two).recycle();

        // recycle the emmigrant
        (*( (*d).getMembers()[x] )).recycle() ;
        // bring in the immigrant
        (*d).immigrate(x, migrant);
    }

}

// SAMPLE MIGRANTS FROM RANDOM DEMES THROUGHOUT THE META-POPULATION -- AS PER KIMURA AND WEISS (1964), M-INFINITY TERM
void MetaPop::longRangeDisp(Deme* d)    // d receives the random deme's emigrant
{
    double lrdRate = 0.01;
    int num_migrants = lrdRate * (*d).getSize(); //gsl_ran_poisson(r, lrdRate);

    for ( int i = 0 ; i < num_migrants ; i++ )
    {
        Deme* e;                    // The randomly chosen deme that supplies the immigrant to d
        Individual* rand_ind_imm;     // the emigrant to d
        Individual* rand_ind_emm;     // the emigrant to e
        int imm;
        int emm;

        do
        {
            e = m_demes[gsl_rng_uniform(s_rand_num) * m_num_demes];
            imm = gsl_rng_uniform(s_rand_num) * (*e).getSize();
            rand_ind_imm = (*e).getMembers()[imm];
        }
        while ( ((*rand_ind_imm).p_migrated) );

        do
        {
            emm = gsl_rng_uniform(s_rand_num) * (*d).getSize();
            rand_ind_emm = (*d).getMembers()[emm];
        }
        while ( ((*rand_ind_emm).p_migrated) );

        Individual* temp_ind = (*d).getInd(emm);
        (*d).immigrate(emm, (*e).getInd(imm) );
        (*e).immigrate(imm, temp_ind);
    }
}

// MAKE A NEW GENERATION OF THE HYBRID ZONE
void MetaPop::newGeneration()
{
    // perform the migrations for the generation
    (*this).stepStone();

    std::vector<Deme*>::iterator iter_d;

    // then run a genertion of mating and reproduction within each deme
    for ( iter_d = m_demes.begin() ; iter_d < m_demes.end() ; iter_d++ )
    {
        (*(*iter_d)).newGeneration();
    }
}

// COLLECT SAMPLES OF SIZE N FROM EACH DEME IN THE METAPOPULATION
Individual** MetaPop::sampleMetapop(int n)
{
    Individual** new_sample = new Individual* [m_num_demes * n];
    std::vector<Deme*>::iterator iter_d;
    int d = 0;

    for( iter_d = m_demes.begin() ; iter_d < m_demes.end() ; iter_d++, d++ )
    {
        for(int i = 0 ; i < n ; i++)
        {
            new_sample[n*d + i] = (**iter_d).randInd();
            (*new_sample[n*d + i]).locateJunctions();
        }
    }

    return new_sample;
}

#endif // METAPOP_H
