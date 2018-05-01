
#ifndef JUNGEN_ACC
#define JUNGEN_ACC

#include <cmath>
#include <vector>
#include <algorithm>
#include <string>
#include <sstream>

#include "junction.h"
#include "gamete.h"
#include "chromosome.h"
#include "autosome.h"
#include "sex_chromosome.h"
#include "cytoplasm.h"

extern std::ofstream g_check_file;


// Take the sum over integers from 1 to n
int sigmaNum(int n)
{
    int sum = 0;
    for ( int i = 1 ; i <= n ; i++ )
    {
        sum += i;
    }
    return sum;
}

// calculate a factorial
int factorial (int n)
{
    int fact = n;
    for( int m = n - 1 ; m > 1 ; m-- )
    {
        fact = fact * m;
    }
    return fact;
}

// partial factorial -- calculate a product over n to k
int part_factorial (int n, int k)
{
    int p_fact = n;
    for ( int m = n - 1 ; m > k ; m-- )
    {
        p_fact = p_fact * m;
    }
    return p_fact;
}

// calculate n choose k
int n_choose_k (int n, int k)
{
    int nck = 0;
    nck = part_factorial(n, n-k) / factorial(k);
    return nck;
}

// express an integer as a string from a different base system -- base < 10
std::vector<int> change_base(int x, int base, int num_digits)
{
    std::vector<int> new_number;

    // find the largest base power to go into the x
    int i = 0;

    do
    {
        ++i;
    } while ( pow(base, i) <= x );

    for( int nd = num_digits ; nd > i ; nd-- )
    {
        new_number.push_back(0);
    }

    --i;

    // count the number of times the next greatest base power fits into the number
    while ( i >= 0 )
    {
        int digit = 0;
        do
        {
            digit++;
        } while ( pow(base, i) * digit <= x );
        new_number.push_back(--digit);

        x = x - ( pow(base, i) * digit );

        --i;
    }

    return new_number;
}

// take a vector of digits representing a number in a base system 'b', and convert it to a decimal number
int convert_to_decimal(std::vector<int> number, int base)
{
    std::vector<int>::iterator iter_d;
    int decimal = 0;
    int i;
    for(iter_d = number.begin(), i = number.size() - 1 ; iter_d < number.end() ; ++iter_d, --i)
    {
        decimal = decimal + (*iter_d) * pow(base, i);
    }
    return decimal;
}

// convert a vector of ints to a string
std::string vec_ints_to_string(std::vector<int> i)
{
    std::stringstream str;
    std::vector<int>::iterator iter_i;
    for( iter_i = i.begin() ; iter_i < i.end() ; iter_i++ )
    {
        str << (*iter_i);
    }

    std::string key = str.str();
    return key;
}



//MAKE A GAMETE POOL FOR A POPULATION WITH XY CHROMOSOMES
std::vector<Gamete*> XYgametePool(int anc, std::vector<double> lengths, std::vector<ChrType> types, double p = 0.001)
{
    std::vector<Gamete*> gamete_pool;

    for ( int i = 0 ; i < anc ; i++ )
    {
        int chr_n = 0;
        std::vector<Chromosome*> chr_pool;
        Chromosome* hetGamChr;
        Junction* cen = new ( Junction::s_newAddyAssign() ) Junction(chr_n, 0, i);
        Junction::s_junction_pool.push_back(cen);
        Junction* tel = new ( Junction::s_newAddyAssign() ) Junction(chr_n, lengths[0], i);
        Junction::s_junction_pool.push_back(tel);
        CNode *cent1 = new ( CNode::newAddyAssign() ) CNode(cen);

        (*cent1).setLeftNode(0);
        CNode *telo1 = new ( CNode::newAddyAssign() ) CNode(tel);
        (*telo1).setRightNode(0);

        (*cent1).setRightNode(telo1);
        (*telo1).setLeftNode(cent1);

        hetGamChr = new SexChromosome(cent1, telo1, 0, chr_n, lengths[0], Y);
        chr_pool.push_back(hetGamChr);


        std::vector<ChrType>::iterator iter_t;
        std::vector<double>::iterator iter_l = lengths.begin();


        for(iter_t = types.begin(); iter_t < types.end() ; iter_t++)
        {
            //errorCheck << "ChrNumber: " << chr_n << " " << *iter_t << endl;

            Junction* cen = new Junction(chr_n, 0, i);
            Junction::s_junction_pool.push_back(cen);
            Junction* tel = new Junction(chr_n, *iter_l, i);
            Junction::s_junction_pool.push_back(tel);

            CNode *cent1 = new (CNode::newAddyAssign() ) CNode(cen);
            (*cent1).setLeftNode(0);
            CNode *telo1 = new (CNode::newAddyAssign() ) CNode(tel);
            (*telo1).setRightNode(0);

            (*cent1).setRightNode(telo1);
            (*telo1).setLeftNode(cent1);

            Chromosome* new_chr;
            if ( *iter_t == X )
            {
                new_chr = new SexChromosome(cent1, telo1, 0, chr_n, *iter_l, *iter_t);
            }
            else if ( (*iter_t) == M || (*iter_t) == C || (*iter_t) == CP )
            {
                new_chr = new Cytoplasm(cent1, telo1, 0, chr_n, 0, *iter_t);
                (*tel).setPosition(0);
            }
            else
            {
                new_chr = new Autosome(cent1, telo1, 0, chr_n, *iter_l, *iter_t);
            }

            //(*new_chr).displayChromosome();

            chr_pool.push_back(new_chr);

            chr_n++;
            iter_l++;
        }

        std::vector<Chromosome*> male_Genome;
        male_Genome.push_back(chr_pool[0]);

        std::vector<Chromosome*>::iterator iter_c;
        for( iter_c = chr_pool.begin() + 2 ; iter_c < chr_pool.end() ; iter_c++ )
        {
            male_Genome.push_back(*iter_c);
        }

        Gamete* male_gamete = new Gamete(male_Genome, 0);
        //(*male_gamete).displayGamete();
        gamete_pool.push_back(male_gamete);

        std::vector<Chromosome*> female_Genome;
        for( iter_c = chr_pool.begin() + 1 ; iter_c < chr_pool.end() ; iter_c++ )
        {
            female_Genome.push_back(*iter_c);
        }

        Gamete* female_gamete = new Gamete(female_Genome, 0);
        //(*female_gamete).displayGamete();
        gamete_pool.push_back(female_gamete);

    }
    return gamete_pool;
}

//MAKE A GAMETE POOL FOR A POPULATION WITH ZW CHROMOSOMES
std::vector<Gamete*> ZWgametePool(int anc, std::vector<double> lengths, std::vector<ChrType> types, double p = 0.001)
{
    std::vector<Gamete*> gamete_pool;

    for ( int i = 0 ; i < anc ; i++ )
    {
        int chr_n = 0;
        std::vector<Chromosome*> chr_pool;
        Chromosome* hetGamChr;

        Junction* cen = new ( Junction::s_newAddyAssign() ) Junction(chr_n, 0, i);
        //(Junction::s_junction_pool).push_back(cen);
        Junction* tel = new ( Junction::s_newAddyAssign() ) Junction(chr_n, lengths[0], i);
        //(Junction::s_junction_pool).push_back(tel);

        CNode *cent1 = new ( CNode::newAddyAssign() ) CNode(cen);
        (*cent1).setLeftNode(0);

        CNode *telo1 = new ( CNode::newAddyAssign() ) CNode(tel);
        (*telo1).setRightNode(0);

        (*cent1).setRightNode(telo1);
        (*telo1).setLeftNode(cent1);

        hetGamChr = new SexChromosome(cent1, telo1, 0, chr_n, lengths[0], W);
        chr_pool.push_back(hetGamChr);


        std::vector<ChrType>::iterator iter_t;
        std::vector<double>::iterator iter_l;
        chr_n++;

        for(iter_t = types.begin(), iter_l = lengths.begin() ; iter_t < types.end() , iter_l < lengths.end() ; iter_t++, iter_l++ )
        {
            cen = new ( Junction::s_newAddyAssign() ) Junction(chr_n, 0, i);
            Junction::s_junction_pool.push_back(cen);
            tel = new ( Junction::s_newAddyAssign() ) Junction(chr_n, *iter_l, i);
            Junction::s_junction_pool.push_back(tel);


            cent1 = new (CNode::newAddyAssign() ) CNode(cen);
            (*cent1).setLeftNode(0);
            telo1 = new (CNode::newAddyAssign() ) CNode(tel);
            (*telo1).setRightNode(0);

            (*cent1).setRightNode(telo1);
            (*telo1).setLeftNode(cent1);

            Chromosome* new_chr;
            if ( *iter_t == Z )
            {
                new_chr = new ( Chromosome::s_newAddyAssign() ) SexChromosome(cent1, telo1, 0, chr_n, *iter_l, *iter_t);
            }
            else if ( (*iter_t) == M || (*iter_t) == C || (*iter_t) == CP )
            {
                new_chr = new ( Chromosome::s_newAddyAssign() ) Cytoplasm(cent1, telo1, 0, chr_n, 0, *iter_t);
            }
            else
            {
                new_chr = new (Chromosome::s_newAddyAssign() ) Autosome(cent1, telo1, 0, chr_n, *iter_l, *iter_t);
            }

            //(*new_chr).displayChromosome();
            Chromosome copyChromosome(*new_chr);

            chr_pool.push_back(new_chr);

            chr_n++;
        }

        std::vector<Chromosome*> female_Genome;
        female_Genome.push_back(chr_pool[0]);

        std::vector<Chromosome*>::iterator iter_c;
        for( iter_c = chr_pool.begin() + 2 ; iter_c < chr_pool.end() ; iter_c++ )
        {
            female_Genome.push_back(*iter_c);
        }


        std::vector<Chromosome*> male_Genome;
        for( iter_c = chr_pool.begin() + 1 ; iter_c < chr_pool.end() ; iter_c++ )
        {
            male_Genome.push_back(*iter_c);
        }

        Gamete* male_gamete = new Gamete(male_Genome, 0);
        gamete_pool.push_back(male_gamete);

        Gamete* female_gamete = new Gamete(female_Genome, 0);
        gamete_pool.push_back(female_gamete);

    }

    return gamete_pool;
}
//MAKE A GAMETE POOL FOR A POPULATION OF DIPLOIDS WITHOUT SEX CHROMOSOMES
// MAKE THE CHROMOSOMES NEEDED TO START THE POPULATION
std::vector<Gamete*> HHgametePool(int anc, std::vector<double> lengths, std::vector<ChrType> types)
{
    std::vector<Gamete*> gamete_pool;

    for ( int i = 0 ; i < anc ; i++ )
    {
        int chr_n = 0;
        std::vector<Chromosome*> chr_pool;

        std::vector<ChrType>::iterator iter_t;
        std::vector<double>::iterator iter_l;
        chr_n++;



        for(iter_t = types.begin(), iter_l = lengths.begin() ; iter_t < types.end() , iter_l < lengths.end() ; iter_t++, iter_l++ )
        {
            Junction* cen = new Junction(chr_n, 0, i);
            Junction* tel = new Junction(chr_n, *iter_l, i);

            CNode *cent1 = new (CNode::newAddyAssign() ) CNode(cen);
            (*cent1).setLeftNode(0);
            CNode *telo1 = new (CNode::newAddyAssign() ) CNode(tel);
            (*telo1).setRightNode(0);

            (*cent1).setRightNode(telo1);
            (*telo1).setLeftNode(cent1);

            Chromosome* new_chr;

            if ( (*iter_t) == M || (*iter_t) == C || (*iter_t) == CP )
            {
                new_chr = new Cytoplasm(cent1, telo1, 0, chr_n, 0, *iter_t);
            }
            else
            {
                new_chr = new Autosome(cent1, telo1, 0, chr_n, *iter_l, *iter_t);
            }

            Chromosome copyChromosome(*new_chr);

            chr_pool.push_back(new_chr);

            chr_n++;
        }

        std::vector<Chromosome*> female_Genome;
        female_Genome.push_back(chr_pool[0]);

        std::vector<Chromosome*>::iterator iter_c;
        for( iter_c = chr_pool.begin() + 2 ; iter_c < chr_pool.end() ; iter_c++ )
        {
            female_Genome.push_back(*iter_c);
        }


        std::vector<Chromosome*> male_Genome;
        for( iter_c = chr_pool.begin() + 1 ; iter_c < chr_pool.end() ; iter_c++ )
        {
            male_Genome.push_back(*iter_c);
        }

        Gamete* male_gamete = new Gamete(male_Genome, 0);
        gamete_pool.push_back(male_gamete);

        Gamete* female_gamete = new Gamete(female_Genome, 0);
        gamete_pool.push_back(female_gamete);

    }

    return gamete_pool;
}

#endif //JUNGEN_ACC
