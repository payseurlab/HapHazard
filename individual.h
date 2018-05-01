
#ifndef INDIVIUDUAL003_H
#define INDIVIUDUAL003_H

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <deque>

#include "chromosome.h"
#include "autosome.h"
#include "sex_chromosome.h"
#include "cytoplasm.h"
#include "junction.h"
#include "gene.h"
#include "landscape.h"
#include "jungen_acc.cpp"

extern std::ofstream g_check_file;

enum MatingSys { HH = 0, XY = 1, ZW = 2 };

/*
    INDIVIDUAL

    Individuals are the most complex objects in the simulator.
    They cause or respond to many factors: genetics, selection, phenotypes, recombination, sampling/summarization... etc.

    Most importantly, genetics and recombination of the genome takes place inside an individual.

*/


class Individual
{
    public:
        // constuctor prototype
        Individual(Gamete* mG, Gamete* pG, int l, bool m, int nP );
        friend class Gene;      // genes need access to several of individual's variables

        // destructor prototype
        ~Individual();

        // getters
        std::vector<Chromosome*> getGenome() { return m_genome; }
        Chromosome* getChr(int i) { return m_genome[i]; }
        int getSex() { return m_sex; }
        int getLocation() { return m_location; }
        double getRFitness() { return m_repr_fitness; }
        double getDFitness() { return m_dev_fitness; }
        double getEFitness() { return m_env_fitness; }
        double getPhenotype(int i) { return m_phenotypes[i]; }
        int getLifespan() { return m_lifespan; }

        // setters
        void setLocation(int l) { m_location = l; }

        // Biological/Genetic Functions/variables
        void determineSex();                                            // checks for Y or Z chromosomes and determines the sex of the individual
        void setRFitness( double r ) { m_repr_fitness = r; }            // sets the reproductive fitness
        void setDFitness( double d ) { m_dev_fitness = d; }             // sets the developmental fitness
        void setEFitness( double e ) { m_env_fitness = e; }             // sets the environmental fitness
        void setPhenotype( int i, double p ) { m_phenotypes[i] = p; }   // sets the value for a phenotype indexed by i
        Chromosome* gameteChromosome(Chromosome* c, Chromosome* d);     // returns one chromosome or the other
        Chromosome* gameteHetChromosome(Chromosome* c, Chromosome* d );
        Gamete* makeGamete();
        std::vector<double> oneCO(double m);
        std::vector<double> coNoInt(int xo, double l);
        std::vector<double> coGamInt(double l);
        std::vector<double> parCrossovers(std::vector<double> co, double p);
        Chromosome * recombination(Chromosome* c, Chromosome* d, std::vector<double> co_list, ChrType t);
        double phenotype(std::vector<Chromosome*> g, Landscape* l);
        void locateJunctions();

        // Summary Statistics Calculators
        int calcnumJunctOne(int x);
        double calcATL(int x);
        double calcHybridIndex();
        double calcFractionHeterogenic();
        double calcFractionHomogenic(int a);
        std::vector<Block> getBlocks(int c);

        // Data Output/Display
        void displayChromosomes();
        void textChromosomes(int x);

        // Memory Management
        void recycle();
        static void s_drainIndividualReservoir();

        // STATIC FUNCTIONS
        static Individual* newAddyAssign();

        bool p_migrated;  // true/false, did the individual migrate this generation

        // Static variables
        static MatingSys s_gp_type;
        static std::vector<Gene*> s_genes;
        static std::vector<Gamete*> s_ind_gamete_pool;
        static gsl_rng* s_rand_num;
        static double s_min_block_size;

        static std::deque<Individual*> s_individual_reservoir;

        static int s_num_chr;
        static int s_interference_opt;
        static int s_num_of_phenotypes;
        static std::vector<Landscape*> s_landscapes;

    private:
        std::vector<double> m_phenotypes;    // a list of phenotype values
        std::vector<Chromosome*> m_genome;   // a list of chromosomes that make up the genome
        double m_repr_fitness;          // the probability of passing on a gamete successfully
        double m_dev_fitness;           // the probability of surviving to adulthood
        double m_env_fitness;           // the probability of surviving in the environment
        int m_sex;                      // the sex of the individual
        int m_lifespan;                 // the generation the individual lived in
        int m_location;                 // the current location of the individual

};

// INDIVIDUAL CONSTRUCTOR
Individual::Individual(Gamete* mat_gamete, Gamete* pat_gamete, int l, bool m = false, int num_phen = 3 ):
    m_lifespan(g_generation), m_location(l)
{
    p_migrated = false;
    m_phenotypes.reserve(s_num_of_phenotypes);

    int gen_size = 2 * (s_num_chr - 1) + 1; // calculate the number of chromosomes (sex chromosomes, autosomes, + 1 cytoplasmic chromosome)

    m_genome.reserve( gen_size ); // reserve space for the genome

    std::vector<Chromosome*> gam1 = (*mat_gamete).getHapGenome();     // copy the parental gametes
    std::vector<Chromosome*> gam2 = (*pat_gamete).getHapGenome();

    // Determine the sex of the individual based on gametes
    // iterate over both gametes' chromosomes simultaneously and check for sex chromosome types
    std::vector<Chromosome*>::iterator iter_g1;
    std::vector<Chromosome*>::iterator iter_g2;

    m_sex = -1;

    for( iter_g1 = gam1.begin() , iter_g2 = gam2.begin() ; iter_g1 < gam1.end() - 1, iter_g2 < gam2.end() - 1 ; iter_g1++, iter_g2++ )
    {
        ChrType chromo_type_one = (**iter_g1).getType();
        ChrType chromo_type_two = (**iter_g2).getType();

        if ( chromo_type_one == Y || chromo_type_two == Y )
        {
            m_sex = 1;
        }
        else if ( chromo_type_one == W || chromo_type_two == W )
        {
            m_sex = 0;
        }
        else if ( chromo_type_one == X && chromo_type_two == X )
        {
            m_sex = 0;
        }
        else if ( chromo_type_one == Z && chromo_type_two == Z )
        {
            m_sex = 1;
        }

        Chromosome* new_chr_one = (**iter_g1).duplicateChr(&m_location);
        Chromosome* new_chr_two = (**iter_g2).duplicateChr(&m_location);

        m_genome.push_back(new_chr_one);
        m_genome.push_back(new_chr_two);

    }

    ChrType cytoplasm_type = (**iter_g1).getType();

    switch(cytoplasm_type)
    {
        case M:
        {
            Chromosome* new_cyto = (**iter_g1).duplicateChr(&m_location);
            m_genome.push_back(new_cyto);
            break;
        }
        case C:
        {
            Chromosome* new_cyto = (**iter_g1).duplicateChr(&m_location);
            m_genome.push_back(new_cyto);
            break;
        }
        case CP:
        {
            Chromosome* new_cyto = (**iter_g2).duplicateChr(&m_location);
            m_genome.push_back(new_cyto);
            break;
        }
	default:
	{
	    std::cout << "Non-cytoplamic chromosome where cytoplasm required.\n" << std::endl;
	    exit(1);
	    break;
	}

    }

    // Determine the phenotypes of the individual
    // for each gene, get the genotype and alter the phenotype accordingly
    // calculate the reproductive fitness
    m_phenotypes[0] = phenotype(m_genome, s_landscapes[0]);

    // NOTE: DEV and ENV fitness calculators commented out because they are not developed
    // but will be used later

    // calculate the developmental fitness
    // double m_dev_fitness = phenotype(m_genome, dFitLand);
    // phenotypes.push_back(m_dev_fitness);

    // calculate the environmental fitness
    // double m_env_fitness = phenotype(m_genome, eFitLand);
    // phenotypes.push_back(m_env_fitness);

}

// INDIVIDUAL DESTRUCTOR
Individual::~Individual()
{
    std::vector<Chromosome*>::iterator iter_c;
    //std::std::cout << "Ind: " << this << std::endl;
    for ( iter_c = m_genome.begin() ; iter_c < m_genome.end() ; iter_c++ )
    {
        (**iter_c).recycle();
    }
}

// INDIVIDUAL RECYCLER
void Individual::recycle()
{
    //std::std::cout << "IndRec" << std::endl;
    p_migrated = false;
    s_individual_reservoir.push_back(this);
    std::vector<Chromosome*>::iterator iter_c;
    //std::std::cout << "Ind: " << this << std::endl;
    for ( iter_c = m_genome.begin() ; iter_c < m_genome.end() ; iter_c++ )
    {
        (**iter_c).recycle();
    }
}

// CLEAR ALL INDIVIDUALS FROM MEMORY
void Individual::s_drainIndividualReservoir()
{
    std::deque<Individual*>::iterator iter_i;
    for(iter_i = s_individual_reservoir.begin() ; iter_i != s_individual_reservoir.end() ; iter_i++ )
    {
        delete *iter_i;
    }
}

// GET AN INDIVIDUAL BLOCK FROM THE INDIVIDUAL RESERVOIR
Individual* Individual::newAddyAssign()
{
    Individual* new_addy;
    if ( s_individual_reservoir.size() != 0 )
    {
        new_addy = *s_individual_reservoir.begin();
        s_individual_reservoir.pop_front();
    }
    else
    {
        new_addy = new Individual(s_ind_gamete_pool[0], s_ind_gamete_pool[0], 0);
    }

    return new_addy;
}

// RETURN A VECTOR OF THE ANCESTRY BLOCKS ON THE INDIVIDUAL'S CHROMOSOMES
std::vector<Block> Individual::getBlocks(int c)
{
    Chromosome* chr_one = m_genome[2*c];
    Chromosome* chr_two = m_genome[2*c+1];

    std::vector<Block> blocks_one = (*chr_one).getBlocks(c);
    std::vector<Block> blocks_two = (*chr_two).getBlocks(c);
    blocks_one.insert(blocks_one.end(), blocks_two.begin(), blocks_two.end() );

    return blocks_one;
}

// UPDATE THE LOCATION COUNTS FOR EACH JUNCTION IN THIS INDIVIDUAL
void Individual::locateJunctions()
{
    std::vector<Chromosome*>::iterator iter_c;
    for( iter_c = m_genome.begin() ; iter_c < m_genome.end() ; iter_c++ )
    {
        (**iter_c).locateJunctions();
    }
}

// PRINT THE CHROMOSOMES TO THE TERMINAL
void Individual::displayChromosomes()
{
    //std::std::cout << this << std::endl;
    int x = m_genome.size();
    for ( int i = 0 ; i < x ; i++ )
    {
        std::cout << m_genome[i] << std::endl;
        (*m_genome[i]).displayChromosome();
    }
}

// PRINT THE CHROMOSOMES TO A TEXT FILE
void Individual::textChromosomes(int x = 2)
{
    std::vector<Chromosome*>::iterator iter_c;
    for ( iter_c = m_genome.begin() ; iter_c < m_genome.end() ; iter_c++ )
    {
        (**iter_c).textChromosome();
    }
    g_check_file << std::endl;
}

// MAKE A GAMETE FROM AN INDIVIDUAL
Gamete* Individual::makeGamete()
{
    std::vector<Chromosome*> gamete;
    Gamete* new_gamete;

    Chromosome* g_chr;

    ChrType c_type1;
    ChrType c_type2;

    for ( unsigned int i = 0 ; i < m_genome.size() ; i++ )
    {

        if ( i % 2 != 0 )
        {
            i++;
        }

        if ( i == m_genome.size() - 1 )
        {
            Chromosome* new_cytoplasm = (*m_genome[i]).duplicateChr(&m_location);
            gamete.push_back( new_cytoplasm );
        }
        else
        {
            c_type1 = (*m_genome[i]).getType();
            c_type2 = (*m_genome[i+1]).getType();


            if( ( c_type1 == A && c_type2 == A ) || ( c_type1 == X && c_type2 == X ) || ( c_type1 == Z && c_type2 == Z ) )
            {
                g_chr = gameteChromosome(m_genome[i], m_genome[i+1] );
                gamete.push_back(g_chr);
            }
            else if ( ( c_type1 == X && c_type2 == Y ) || ( c_type1 == Y && c_type2 == X ) || ( c_type1 == W && c_type2 == Z ) || ( c_type1 == Z && c_type2 == W )   )
            {
                g_chr = gameteHetChromosome( m_genome[i], m_genome[i+1] );
                gamete.push_back(g_chr);
            }
            else
            {
                std::cout << "Mis-Ordered genomic vector. Must be S,S,A1,A1,A2,A2,...,An,An,C" << std::endl;
                exit(1);
            }
        }
    }

    new_gamete = new Gamete(gamete, &m_location);

    return new_gamete;
}

// MAKE A CHROMOSOME FOR THE GAMETE
Chromosome* Individual::gameteChromosome(Chromosome* c, Chromosome* d)
{
    double mu = (*c).getLength();               // the genetic length of the chromosome provides the mean of a Poisson distribution of crossover number
    double strand_opt = gsl_rng_uniform(s_rand_num);      // strand_opt chooses which of the parental gametes will be chosen, either to be passed on, or as the starting strand for recombination
    Chromosome* gam_chr;                        // a vector of pointers to the chromosomes that will compose the gamete
    std::vector<double> crossovers;                  // this vector stores a list of crossovers to be used by the recombination function if needed

    // make the gamete according to the type of recombination being used
    // 0 : one crossover per chromosome
    // 1 : poisson number of crossovers per chromosome
    // 2 : gamma interference between crossovers
    switch(s_interference_opt)
    {
        case 0:
        {
            // ONE CROSSOVER PER CHROMOSOME
            // first randomly decide whether to recombine
            int rec_choice = gsl_rng_uniform(s_rand_num) + 0.5;
            if(rec_choice == 0) // if not return the chromosome chosen by strand_opt
            {
                if ( strand_opt < 0.5 )
                {
                    gam_chr = (*c).duplicateChr(&m_location);
                }
                else
                {
                    gam_chr = (*d).duplicateChr(&m_location);
                }
            }
            else    // otherwise return a recombinant, whose left end comes from the chromosome chosen by strand_opt
            {
                double one_crossover = gsl_rng_uniform(s_rand_num) * mu; // pick a crossover
                crossovers.push_back(one_crossover); // add it to the crossover list
                if ( strand_opt < 0.5 )
                {
                    gam_chr = recombination(c, d, crossovers, (*d).getType() );
                }
                else
                {
                    gam_chr = recombination(d, c, crossovers, (*c).getType() );
                }
            }
            break;
        }
        case 1:
        {
            // POISSON CROSSOVERS
            // first pick the number of crossovers from a poisson distribution
            // with the mean/variance equal to the genetic length of the chromosome
            int x = gsl_ran_poisson(s_rand_num, mu);
            if(x == 0) // if the number of crossovers is 0, do not recombine, just return the strand_opt
            {
                if ( strand_opt < 0.5 )
                {
                    gam_chr = (*c).duplicateChr(&m_location);
                }
                else
                {
                    gam_chr = (*d).duplicateChr(&m_location);
                }
            }
            else // if the crossver number is 1+, recombine
            {
                crossovers = coNoInt(x, mu); // generate the crossovers
                if ( strand_opt < 0.5 )
                {
                    gam_chr = recombination(c, d, crossovers, (*d).getType() );
                }
                else
                {
                    gam_chr = recombination(d, c, crossovers, (*c).getType() );
                }
            }
            break;
        }
        case 2:
        {
            // Gamm interference of crossovers
            int x = gsl_rng_uniform(s_rand_num); // choose whether or not to recombine
            if (x < 0.5)
            {
                if ( strand_opt < 0.5 )
                {
                    gam_chr = (*c).duplicateChr(&m_location);
                }
                else
                {
                    gam_chr = (*d).duplicateChr(&m_location);
                }
            }
            else  // otherwise
            {
               crossovers = coGamInt(mu);   // generate a crossover list with gamma interference
                if ( strand_opt < 0.5 )
                {
                    gam_chr = recombination(c, d, crossovers, (*d).getType() );
                }
                else
                {
                    gam_chr = recombination(d, c, crossovers, (*c).getType() );
                }
            }
            break;
        }
        default:
        {
            std::cout << "Recombination interference option not specified correctly (0, 1, or 2)!" << std::endl;
            exit(1);
            break;
        }
    }

    return gam_chr;
}


// GET A SEX CHROMOSOME WHEN THE SEX CHROMOSOMES ARE HETEROGAMETIC
Chromosome* Individual::gameteHetChromosome(Chromosome* c, Chromosome* d)
{
    double mu;
    double par_boundary;

    if ( (*c).getType() == Y || (*c).getType() == W )
    {
        mu = (*c).getLength();        // the genetic length of the chromosome provides the mean of a Poisson distribution of crossover number
        par_boundary = (*c).getParB();
    }
    else
    {
        mu = (*d).getLength();
        par_boundary = (*d).getParB();
    }

    double strand_opt = gsl_rng_uniform(s_rand_num);      // strand_opt chooses which of the parental gametes will be chosen, either to be passed on, or as the starting strand for recombination
    Chromosome* gam_chr;                        // a vector of pointers to the chromosomes that will compose the gamete
    std::vector<double> crossovers;                  // this vector stores a list of crossovers to be used by the recombination function if needed

    // allow only one crossover on the par
    int x = gsl_rng_uniform(s_rand_num) + 0.5;

    if(x == 0)
    {
        if ( strand_opt < 0.5 )
        {
            gam_chr = (*c).duplicateChr(&m_location);
        }
        else
        {
            gam_chr = (*d).duplicateChr(&m_location);
        }
    }
    else
    {
        double one_crossover = gsl_rng_uniform(s_rand_num) * mu * par_boundary;
        crossovers.push_back(one_crossover);

        if ( strand_opt < 0.5 )
        {
            gam_chr = recombination(c, d, crossovers, (*d).getType() );
        }
        else
        {
            gam_chr = recombination(d, c, crossovers, (*c).getType() );
        }
    }


    return gam_chr;
}

// RECOMBINE THE TWO CHROMOSOMES OF AN INDIVIDUAL
Chromosome* Individual::recombination(Chromosome* c, Chromosome* d, std::vector<double> co_list, ChrType t)
{
    //std::cout << "RR_AA" << std::endl;
    CNode* cur_a = 0;     // cur_a is a pointer to the following CNode on strand that is currently being copied to the recombinant
    CNode* cur_b = 0;     // cur_b is the leading CNode
    CNode* opp_a = 0;     // the follower on the opposite strand
    CNode* opp_b = 0;     // the leader on the opposite strand

    cur_b = (*c).getLeftEnd();
    opp_b = (*d).getLeftEnd();

    CNode* rec_a = new ( (*cur_b).newAddyAssign() ) CNode( (*cur_b).getJunction() );        // rec_a is the following CNode pointer on the recombinant, and starts by pointing to the centromere on the current strand

    CNode* rec_f;

    cur_b = (*cur_b).getRightNode();                              // advance cur_b to the next CNode

    Chromosome* recombinant = 0;
    switch(t)
    {
        case A:
        {
            recombinant = new ( Chromosome::s_newAddyAssign() ) Autosome(rec_a, 0, &m_location, (*c).getNumber(), (*c).getLength(), t );      // initialize a recombinant chromosome, and give it the centromere, rec_a, and leave the telomere as NULL
            break;
        }
        case W:
        {
            recombinant = new ( Chromosome::s_newAddyAssign() ) SexChromosome(rec_a, 0, &m_location, (*c).getNumber(), (*c).getLength(), t );      // initialize a recombinant chromosome, and give it the centromere, rec_a, and leave the telomere as NULL
            break;
        }
        case X:
        {
            recombinant = new ( Chromosome::s_newAddyAssign() ) SexChromosome(rec_a, 0, &m_location, (*c).getNumber(), (*c).getLength(), t );      // initialize a recombinant chromosome, and give it the centromere, rec_a, and leave the telomere as NULL
            break;
        }
        case Y:
        {
            recombinant = new ( Chromosome::s_newAddyAssign() ) SexChromosome(rec_a, 0, &m_location, (*c).getNumber(), (*c).getLength(), t );      // initialize a recombinant chromosome, and give it the centromere, rec_a, and leave the telomere as NULL
            break;
        }
        case Z:
        {
            recombinant = new ( Chromosome::s_newAddyAssign() ) SexChromosome(rec_a, 0, &m_location, (*c).getNumber(), (*c).getLength(), t );      // initialize a recombinant chromosome, and give it the centromere, rec_a, and leave the telomere as NULL
            break;
        }
        default:
        {
            recombinant = new ( Chromosome::s_newAddyAssign() ) Autosome(rec_a, 0, &m_location, (*c).getNumber(), (*c).getLength(), t);      // initialize a recombinant chromosome, and give it the centromere, rec_a, and leave the telomere as NULL
            break;
        }
    }

    std::vector<double>::iterator iter_co;

    for( iter_co = co_list.begin() ; iter_co < co_list.end() ; iter_co++ )   // for each crossover in the list
    {

        while( ( (*cur_b).getJPosition() <= *iter_co && (*cur_b).getRightNode() != 0 ) )  // until cur_b goes past the crossover, or reaches the end of the chromosome
        {
            CNode * rec_b = new ( CNode::newAddyAssign() ) CNode( (*cur_b).getJunction() );    // make a new CNode for the recombinant, and point it to the junction of cur_b's CNode
            (*rec_a).setRightNode(rec_b);                               // make rec_a point to rec_b distally
            (*rec_b).setLeftNode(rec_a);                               // make rec_b point to rec_a proximally
            cur_b = (*cur_b).getRightNode();                            // advance cur_b
            rec_a = (*rec_a).getRightNode();                            // advance rec_a
        }

        cur_a = (*cur_b).getLeftNode();                                // advance cur_a to the node that preceeds cur_b

        while( (*opp_b).getJPosition() <= *iter_co && (*opp_b).getRightNode() != 0 )   // move opp_b past the crossover but not off the chromosome
        {
            opp_b = (*opp_b).getRightNode();
        }

        opp_a = (*opp_b).getLeftNode();              // advance opp_b to the node preceeding opp_b

        int cur_anc = (*cur_a).getJAncestry();    // get the ancestries of the strands at the crossover's position

        int opp_anc = (*opp_a).getJAncestry();


        if ( cur_anc != opp_anc )                 // if the ancestries differ make a new junction and add it to the strand
        {
            double for_block_size = (*opp_b).getJPosition() - (*iter_co);
            double rev_block_size = (*iter_co) - (*cur_a).getJPosition() ;
            //g_check_file << cur_anc << " " << opp_anc << " " << s_min_block_size << " " << for_block_size << " " << rev_block_size << " " ;

            if( for_block_size < s_min_block_size && rev_block_size < s_min_block_size )
            {
                //g_check_file << "BB " << std::endl;
                // if both of the blocks flanking the new junction are too small, add the new junction, but eliminate the two smaller flanking blocks
                Junction * junct = new ( Junction::s_newAddyAssign() ) Junction( (*c).getNumber(), *iter_co, opp_anc); // make a new junction at the crossovers position, with the opposite strand's ancestry

                Junction::s_junction_pool.push_back(junct);        // add then new junction to the junction pool
                CNode *rec_b = new ( CNode::newAddyAssign() ) CNode(junct );                         // make a new CNode that points to the new junction

                CNode* recC;

                if( (*rec_a).getLeftNode() != 0 )   // only get the previous node for rec_a if rec_a is not the beginning of the chromosome
                {
                    recC = rec_a;
                    rec_a = (*rec_a).getLeftNode();
                    (*rec_a).setRightNode(rec_b);
                    (*rec_b).setLeftNode(rec_a);
                    (*recC).recycle();
                }
                else
                {
                    (*rec_a).setRightNode(rec_b);
                    (*rec_b).setLeftNode(rec_a);
                }                                         // recycle the rec_a node

                if( (*opp_b).getRightNode() != 0 )
                {
                    // Only advance if the next block is not the end of the chromosome
                    opp_b = (*opp_b).getRightNode();
                }                                                                     // Advance opp_b one more node before connecting it to the new junction
                rec_b = new ( CNode::newAddyAssign() ) CNode( (*opp_b).getJunction() );                // make a new node to the next junction which is now on the opposite strand
                rec_a = (*rec_a).getRightNode();                                                                      // advance rec_a
                (*rec_a).setRightNode(rec_b);                                                                         // link rec_a and rec_b as before
                (*rec_b).setLeftNode(rec_a);
                rec_f = rec_b;

            }
            else if ( rev_block_size < s_min_block_size && for_block_size > s_min_block_size )
            {
                // if the block behind the new junction is too small, get rid of it without adding a junction
                //g_check_file << "R ";
                CNode *rec_b = new ( CNode::newAddyAssign() ) CNode( (*opp_b).getJunction() );       // skip making and linking a new junciton and crosslink the strands about the new crossover
                //g_check_file << "I_REC_675\t" << *iter_co << std::endl;
                CNode* recC;

                if( (*rec_a).getLeftNode() != 0 )
                {
                    recC = rec_a;
                    rec_a = (*rec_a).getLeftNode();
                    (*rec_a).setRightNode(rec_b);
                    (*rec_b).setLeftNode(rec_a);
                    (*recC).recycle();
                }
                else
                {
                    (*rec_a).setRightNode(rec_b);
                    (*rec_b).setLeftNode(rec_a);
                }

                rec_f = rec_b;

            }
            else if ( for_block_size < s_min_block_size && rev_block_size > s_min_block_size )
            {
                //g_check_file << "F " ;
                // if the block in front of the new junction is too small, get rid of it without adding a junction
                //g_check_file << "I_REC_699\t" << *iter_co << std::endl;
                if( (*opp_b).getRightNode() != 0 )
                {
                    // Only advance if the next block is not the end of the chromosome
                    opp_b = (*opp_b).getRightNode();
                }

                CNode *rec_b = new ( CNode::newAddyAssign() ) CNode( (*opp_b).getJunction() );       // skip making and linking a new junciton and crosslink the strands about the new crossover
                (*rec_a).setRightNode(rec_b);
                (*rec_b).setLeftNode(rec_a);
                rec_f = rec_b;
            }
            else
            {
                //g_check_file << "J ";
                Junction * junct = new ( Junction::s_newAddyAssign() ) Junction( (*c).getNumber(), *iter_co, opp_anc); // make a new junction at the crossovers position, with the opposite strand's ancestry
                //g_check_file << "I_REC_714\t" << *iter_co << std::endl;
                //g_check_file << "NewAddyRec717" << std::endl;
                //(*junct).textJunction();
                Junction::s_junction_pool.push_back(junct);                          // add then new junction to the junction pool
                CNode *rec_b = new ( CNode::newAddyAssign() ) CNode(junct ); // make a new CNode that points to the new junction
                //g_check_file << "NewCNodeRec721" << std::endl;
                //(*junct).textJunction();
                (*rec_a).setRightNode(rec_b);                                 // link rec_a a to the new Junction distally
                (*rec_b).setLeftNode(rec_a);                                 // link the new junction to rec_a proximally
                rec_b = new ( CNode::newAddyAssign() ) CNode( (*opp_b).getJunction() );              // make a new node to the next junction which is now on the opposite strand
                rec_a = (*rec_a).getRightNode();                              // advance rec_a
                (*rec_a).setRightNode(rec_b);                                 // link rec_a and rec_b as before
                (*rec_b).setLeftNode(rec_a);
                rec_f = rec_b;
            }
            //g_check_file << std::endl;
        }
        else                                                        // otherwise
        {
            //std::std::cout << "N";
            CNode *rec_b = new ( CNode::newAddyAssign() ) CNode( (*opp_b).getJunction() );       // skip making and linking a new junciton and crosslink the strands about the new crossover
            (*rec_a).setRightNode(rec_b);
            (*rec_b).setLeftNode(rec_a);
            rec_f = rec_b;
        }
        //std::std::cout << std::endl;
        std::swap(cur_a,opp_a);            // swap the pointers to reflect the new strand orientation when past the crossover
        std::swap(cur_b,opp_b);

    }

    while( (*cur_b).getRightNode() != 0 )                        // after the last crossover, finish the strand up until the telomere
    {
        cur_b = (*cur_b).getRightNode();
        CNode * rec_b = new ( CNode::newAddyAssign() ) CNode( (*cur_b).getJunction() );
        rec_a = (*rec_a).getRightNode();
        (*rec_a).setRightNode(rec_b);
        (*rec_b).setLeftNode(rec_a);
        rec_f = rec_b;
    }

    (*recombinant).setRightEnd(rec_f);
    return recombinant;
}

// ONE CROSSOVER -- GET ONE RANDOM CROSSOVER
std::vector<double> Individual::oneCO(double m)
{
    std::vector<double> co;

    double x = gsl_rng_uniform(s_rand_num) * m;
    co.push_back(x);

    return co;
}

// POISSON CROSSOVERS -- MAKE A LIST OF CROSSOVERS WITHOUT INTERFERENCE
std::vector<double> Individual::coNoInt(int c, double m)
{
    std::vector<double> co;

    for ( int i = 0 ; i < c ; i++)
    {
        double x = gsl_rng_uniform_pos(s_rand_num) * m;
        co.push_back(x);
    }
    sort(co.begin(), co.end());
    return co;
}

// GAMMA CROSSOVERS -- MAKE A LIST OF CROSSOVERS WITH INTERFERENCE
std::vector<double> Individual::coGamInt(double m)
{
    std::vector<double> crossover_list;
    double crossover = m - ( gsl_rng_uniform_pos(s_rand_num) * m );  // the first position
    crossover_list.push_back(crossover);
    double previous_crossover = crossover;
    double dir_choice = 1 - ( gsl_rng_uniform_pos(s_rand_num) );      // choose the direction with the greatest oppurtunity for crossing over

    double v = 8;               // The interference parameter, 8 is the mouse genomewide average as determined in Broman et al. 2002.
    double w = 1 / ( 2 * v );    // the rate parameter for the gamma model from above

    do                                                       // put new crossovers that fall on the chromosome onto the list of crossovers
    {
        double gamma_rand = gsl_ran_gamma(s_rand_num, v, w);              // add crossover chosen from gamma distribution  Function: double gsl_ran_gamma (const gsl_rng * r, double a, double b)
        if (dir_choice < 0.5)
        {
            crossover =  gamma_rand + previous_crossover;
        }
        else
        {
            crossover =  previous_crossover - gamma_rand;
        }

        if ( crossover > 0 && crossover < m )                // if the new crossover is on the chromosome
        {
            crossover_list.push_back(crossover);              // add it to the list
        }
        previous_crossover = crossover;                       // set the new crossover as the previous
    }  while ( crossover > 0 && crossover < m );               // if the new one is off the chromosome, quit adding crossovers

    sort(crossover_list.begin(), crossover_list.end());                // sort the array of positions

    return crossover_list;
}

// REDUCE CROSSOVER POSITIONS TO FIT INTO PSEUDOAUTOSOMAL REGION
std::vector<double> Individual::parCrossovers(std::vector<double> co, double p)
{
    std::vector<double>::iterator iter_c;
    for ( iter_c = co.begin() ; iter_c < co.end() ; iter_c++ )
    {
        *iter_c = *iter_c * p;
    }
    return co;
}

// SUMMARY STATISTICS CALCULATORS

// CALCULATE THE NUMBER OF JUNCTIONS ON ONE CHROMOSOME INDICATED BY X -- 1 = chrOne, 2 = chrTwo, 3 = both chrs
int Individual::calcnumJunctOne(int x)
{
    int num_junct = (*m_genome[x]).calcNumJunctions();
    return num_junct;
}

// CALCULATE THE AVERAGE TRACT LENGTH FOR THE GENOME
double Individual::calcATL(int x)
{
    double ATL = (*m_genome[x]).calcATL();
    return ATL;
}

// CALCULATE THE HYBRID INDEX OF THE GENOME
double Individual::calcHybridIndex()
{
    double hyb_index = -1;

    return hyb_index;
}

// GET THE GENOTYPE KEY FOR A GENE IN AN INDIVIDUAL
double Individual::phenotype(std::vector<Chromosome*> g, Landscape* l)
{

    std::vector<int> genotype;
    std::vector<Gene*> genes = (*l).getLoci();
    std::vector<Gene*>::iterator iter_g;

    for( iter_g = genes.begin() ; iter_g < genes.end() ; iter_g++ )
    {

        unsigned int chr =  (**iter_g).getChr() * 2;
        double pos =  (**iter_g).getPos();

        int anc_one = (*g[ chr     ] ).positionAnc(pos);
        genotype.push_back( anc_one );
        int anc_two = -1;
        if(chr == m_genome.size() - 1 )
        {
            anc_two = anc_one;
        }
        else
        {
            anc_two = (*g[ chr + 1 ] ).positionAnc(pos);
        }

        genotype.push_back( anc_two );

    }
    std::string str_key = vec_ints_to_string(genotype);
    double phenotype = (*l).findPhenotype(str_key);

    return phenotype;
}

#endif // INDIVIDUAL003_H
