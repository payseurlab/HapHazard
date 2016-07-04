#ifndef GAMETE_H
#define GAMETE_H

#include <deque>
#include <vector>

#include "chromosome.h"


class Gamete
{
    public:
        //Gamete constructor prototype
        Gamete(std::vector<Chromosome*> h, int* l);
        // Gamete destructor protoype
        ~Gamete();

        // Make a copy of the gamete
        Gamete* copyGamete(int* l);

        //Getters
        Chromosome* getChromosome(int i) { return m_hap_genome[i]; }
        std::vector<Chromosome*> getHapGenome() { return m_hap_genome; }
        //Chromosome* getCytoplasm() { return cyto; }

        //Setters

        // Data display
        void displayGamete();

        // STATIC FUNCTIONS AND VARIABLES
        // memory management

        void recycle();
        void drainGameteRes();

        static int s_num_gametes;
        static std::deque<Gamete*> s_gamete_reservoir;

    protected:
        std::vector<Chromosome*> m_hap_genome;
		int* m_location;


};

// GAMETE CONSTRUCTOR
Gamete::Gamete(std::vector<Chromosome*> h, int* l):
    m_hap_genome(h), m_location(l)
{
    s_num_gametes++;
}

// MAKE A COPY OF THE GAMETE
Gamete* Gamete::copyGamete(int* l)
{
    std::vector<Chromosome*>::iterator iter_c;
    std::vector<Chromosome*> new_hap_genome;
    for( iter_c = m_hap_genome.begin() ; iter_c != m_hap_genome.end() ; iter_c++ )
    {
        Chromosome* new_chr = (**iter_c).duplicateChr(l);
        new_hap_genome.push_back(new_chr);
    }

    Gamete* new_gamete = new Gamete(new_hap_genome, l);

    return new_gamete;
}

// GAMETE DESTRUCTOR
Gamete::~Gamete()
{
    std::vector<Chromosome*>::iterator iter_c;
    for(iter_c = m_hap_genome.begin() ; iter_c < m_hap_genome.end() ; iter_c++ )
    {
        (**iter_c).recycle();
    }
    s_num_gametes--;
}

// GAMETE RECYCLER
void Gamete::recycle()
{
    std::vector<Chromosome*>::iterator iter_c;
    for(iter_c = m_hap_genome.begin() ; iter_c < m_hap_genome.end() ; iter_c++ )
    {
        (**iter_c).recycle();
    }

    s_gamete_reservoir.push_back(this);
    s_num_gametes--;
}

// CLEAR RESERVED GAMETES FROM MEMORY
void Gamete::drainGameteRes()
{
    std::deque<Gamete*>::iterator iter_g;
    for( iter_g = s_gamete_reservoir.begin() ; iter_g != s_gamete_reservoir.end() ; iter_g++ )
    {
        delete *iter_g;
    }
}

// DISPLAY A GAMETE IN THE TERMINAL
void Gamete::displayGamete()
{
    std::vector<Chromosome*>::iterator iter_g;
    std::cout << "Gamete: " << this << std::endl;

    for( iter_g = m_hap_genome.begin() ; iter_g < m_hap_genome.end() ; iter_g++ )
    {
        std::cout << *iter_g << std::endl;
        (**iter_g).displayChromosome();
    }

}

#endif // GAMETE_H
