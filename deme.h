
#ifndef DEME_H
#define DEME_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <ctime>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "individual.h"
#include "chromosome.h"
#include "jungen_acc.cpp"

extern std::ofstream g_check_file;

/*
    DEME

    Demes represent populations, containing individuals among which random mating can occur.
    They are the sub-populations within a metapopulation, and migration of individuals occurs
    between demes.

    The deme class in this program conatains individuals, and regulates the manner of reproduction.
    The number of individuals is constant through generations
*/

class Deme
{
    public:
        // constructor prototype
        Deme(int s, std::vector<double> a, int l, std::vector<Gamete*> gP);
        ~Deme();

        // getters
        int getSize() { return m_size; }
        std::vector<double> getAncProp() { return m_anc_prop; }
        std::vector<Individual*> getMembers() { return m_members; }
        int* getLocation() { return &m_deme_location;}
        Individual* getInd(int i) { return m_members[i]; }

        // setters
        void immigrate(int i, Individual* j);

        // population functions
        void newGeneration();
        void generationSummary();
        Individual* randInd();
        Individual* randMale();
        Individual* randFemale();
        void purgeMembers();

        // summary statistic calculators
        std::vector<Individual*> sampleDeme(int n);

        // data management
        void textMembers(const char* ofname);

        // static variables
        static gsl_rng* s_rand_num;
        static MatingSys s_mat_sys;
        static int s_num_anc;

    private:
        int m_size;
        std::vector<double> m_anc_prop;
        int m_deme_location;
        std::vector<Individual*> m_members;
        int m_gen;
};

// DEME CONSTRUCTOR
Deme::Deme(int s, std::vector<double> a, int l, std::vector<Gamete*> gP):
    m_size(s), m_anc_prop(a), m_deme_location(l)
{
    // Upon creation in generation 0, we will generate pure individuals according
    // to the propotion of ancestry given... the ancestries are chosen randomly from
    // a multinomial distribution according to the proportions
    unsigned int n[s_num_anc];

    gsl_ran_multinomial(s_rand_num, s_num_anc, m_size, &m_anc_prop[0], n);
    m_gen = g_generation;

    for( int i = 0 ; i < s_num_anc ; i++ )
    {
        for(unsigned int j = 0 ; j <= n[i] ; j++ )
        {
            // Depending on the mating system required, we choose gametes from the pure populations
            // and creat new individuals
            if (Deme::s_mat_sys == XY || Deme::s_mat_sys == ZW )
            {
                //Individual(Gamete* matGamete, Gamete* patGamete, int l = 0, bool m = false ):
                int sex_chr = gsl_ran_bernoulli(s_rand_num, 0.5);
                int G1 = 2 * i + 1;
                int G2 = 2 * i + 1 - sex_chr;

                Gamete* new_mat_gam = (*(Individual::s_ind_gamete_pool[G1])).copyGamete(&m_deme_location);

                Gamete* new_pat_gam = (*(Individual::s_ind_gamete_pool[G2])).copyGamete(&m_deme_location);

                Individual* new_ind = new ( Individual::newAddyAssign() ) Individual( new_mat_gam, new_pat_gam, m_deme_location );

                m_members.push_back(new_ind);

                (*new_mat_gam).recycle();
                (*new_pat_gam).recycle();

            }
            else
            {
                //Individual(Gamete* matGamete, Gamete* patGamete, int l = 0, bool m = false ):
                int G1 = i;
                int G2 = i;
                Gamete* new_mat_gam = (*(Individual::s_ind_gamete_pool[G1])).copyGamete(&m_deme_location);
                Gamete* new_pat_gam = (*(Individual::s_ind_gamete_pool[G2])).copyGamete(&m_deme_location);

                Individual* new_ind = new ( Individual::newAddyAssign() ) Individual( new_mat_gam, new_pat_gam, m_deme_location );

                m_members.push_back(new_ind);

                (*new_mat_gam).recycle();
                (*new_pat_gam).recycle();
            }
        }
    }

}

//DEME DESTRUCTOR
Deme::~Deme()
{
    // when the deme is destroyed, every individual is recycled
    std::vector<Individual*>::iterator iter_i;
    for ( iter_i = m_members.begin() ; iter_i < m_members.end() ; iter_i++ )
    {
        (**iter_i).recycle();
    }
}

// MAKE A NEW GENERATION IN A DEME FROM THE PREVIOUS GENERATION
void Deme::newGeneration()
{
    // each generation results in replacement of all old members with the new
    std::vector<Individual*> newMembers;

    // count through to the populaiton carrying capacity
    for ( int i = 0; i < m_size ; i++ )
    {
        // parents of each new member will be selected randomly
        Individual* parent_one = 0;
        Individual* parent_two = 0;
        bool par_one_pass = false;  // the parents must be tested against each other to ensure
        bool par_two_pass = false;  // male-female pairs, and can be rejected based on fitness

        // select parents at random, then accept or reject them in accord with their reproductive fitness
        do
        {
            double reject_one = gsl_rng_uniform(s_rand_num);

            switch(Deme::s_mat_sys)
            {
                case(HH):
                    parent_one = randInd();
                    break;
                case(XY):
                    parent_one = randFemale();
                    break;
                case(ZW):
                    parent_one = randMale();
                    break;
            }

            //std::cout << "DM_162: " << (*parent_one).getPhenotype(0) << "\t" << reject_one << std::endl;
            if ( reject_one < (*parent_one).getPhenotype(0) )
            {
                par_one_pass = true;
            }

        } while ( par_one_pass == false );

        do
        {

            double reject_two = gsl_rng_uniform(s_rand_num);
            switch(Deme::s_mat_sys)
            {
                case(HH):
                    parent_two = randInd();
                    break;
                case(XY):
                    parent_two = randMale();
                    break;
                case(ZW):
                    parent_two = randFemale();
                    break;
            }


            if ( reject_two < (*parent_two).getPhenotype(0) )
            {
                par_two_pass = true;
            }

        } while ( par_two_pass == false );

        // once parents have been selected, make a gamete from each, and use them to make a new individual
        Gamete* gamete_one = (*parent_one).makeGamete();
        Gamete* gamete_two = (*parent_two).makeGamete();

        Individual* offspring = new ( Individual::newAddyAssign() ) Individual(gamete_one, gamete_two, m_deme_location);

        // Recycle the gametes
        (*gamete_one).recycle();
        (*gamete_two).recycle();

        // Add the new member to the deme
        newMembers.push_back(offspring);

    }

    // once the new generation is complete, delete the old
    purgeMembers();
    m_members = newMembers;
}

// IMMIGRATE A NEW INDIVIDUAL -- add a new indivual from elsewhere
void Deme::immigrate(int i, Individual* j)
{
    m_members[i] = j;
    (*m_members[i]).setLocation(m_deme_location);
    (*m_members[i]).p_migrated = true;              // set migrated to true so it won't migrate again this generation
}

// CHOOSE AN INDIVIDUAL AT RANDOM FROM THE DEME
Individual* Deme::randInd()
{
    int x = gsl_rng_uniform(s_rand_num) * m_size;
    return m_members[x];
}

// SAMPLE N RANDOM INDIVIDUALS FORM THE DEME -- WITHOUT REPLACEMENT
std::vector<Individual*> Deme::sampleDeme( int n)
{
    std::vector<Individual*> sample;

    // sample random individuals from the population until a sample of N is achieved
    for( int i = 0; i < n ; i++ )
    {
        Individual* s_Ind = randInd();

        // before individuals are added to the list, make sure they have not already been sampled
        bool resample = false;
        do
        {
            std::vector<Individual*>::iterator iter_s;
            for( iter_s = sample.begin() ; iter_s < sample.end() ; iter_s++ )
            {
                if( s_Ind == (*iter_s) )
                {
                    resample = true;
                }
            }
        } while (resample == true);

        sample.push_back(s_Ind);
    }
    return sample;
}

// CHOOSE AN INDIVIDUAL AT RANDOM FROM THE DEME
Individual* Deme::randMale()
{
    int x;
    do
    {
        x = gsl_rng_uniform(s_rand_num) * m_size;
    } while ( (*m_members[x]).getSex() != 1 );

    return m_members[x];
}

// CHOOSE AN INDIVIDUAL AT RANDOM FROM THE DEME
Individual* Deme::randFemale()
{
    int x;
    do
    {
        x = gsl_rng_uniform(s_rand_num) * m_size;
    } while ( (*m_members[x]).getSex() != 0 );

    return m_members[x];
}

// PURGE THE OLD DEME FROM MEMORY
void Deme::purgeMembers()
{
    std::vector<Individual*>::iterator iter_m;
    for( iter_m = m_members.begin() ; iter_m < m_members.end() ; iter_m++ )
    {
        (**iter_m).recycle();
    }
}

// PRINT THE MEMORY ADDRESSES OF DEME MEMBERS TO A TEXT FILE
void Deme::textMembers(const char* ofname)
{
    std::fstream deme_out;
    deme_out.open(ofname);

    std::vector<Individual*>::iterator iter_m;
    for( iter_m = m_members.begin() ; iter_m < m_members.end() ; iter_m++ )
    {
        deme_out << (*iter_m) << std::endl;
    }

    deme_out << "*****************" << std::endl;
}

#endif // DEME_H
