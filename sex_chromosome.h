
#ifndef SEXCHROMOSOME_H
#define SEXCHROMOSOME_H

#include "chromosome.h"

extern std::ofstream g_check_file;
extern int g_generation;

class SexChromosome : public Chromosome
{
    public:
        // constructor prototype
        SexChromosome(CNode * left, CNode * right, int* loc, int n, double l, ChrType t);
        // copy constructor
        SexChromosome(SexChromosome& s);
	// destructor
	virtual ~SexChromosome();

        // getters
        double virtual getParB() { return s_par_boundary; }
        double virtual getLength() { return m_length;}
        ChrType virtual getType() { return m_type; }
        int virtual getNumber() { return m_number; }

        //Sex chromosome functions
        virtual int positionAnc(double p);

        // display functions
        void virtual displayChromosome();
        void virtual textChromosome();

        virtual Chromosome* duplicateChr(int* l);

    private:

};

// SEX CHROMOSOME CONSTRUCTOR
SexChromosome::SexChromosome(CNode * cen, CNode * tel, int* loc, int n, double l, ChrType t = X):
    Chromosome(cen, tel, loc, n, l, t)
{

}

// CHROMOSOME COPY CONSTRUCTOR
SexChromosome::SexChromosome(SexChromosome& a):
    Chromosome(a)
{

}

// SEX CHROMOSOME DESTRUCTOR
SexChromosome::~SexChromosome()
{

}

// DUPLICATE AN SEX CHROMOSOME TO BE INHERITED
Chromosome* SexChromosome::duplicateChr(int* l)
{
    CNode * new_left = new  ( CNode::newAddyAssign() ) CNode( (*m_left_end).getJunction() );          // make a new CNode that is a copy of the centromere
    CNode * new_right = new ( CNode::newAddyAssign() ) CNode( (*m_right_end).getJunction() );            // copy the telomere
    Chromosome * rep_chr = new (Chromosome::s_newAddyAssign() ) SexChromosome(new_left, new_right, l, m_number, m_length, m_type);       // make a new chromosome using the copied centromere and telomere

    CNode * temp_node = (*m_left_end).getRightNode();                        // make pointer that starts on the CNode after the centromere on the template chromosome
    CNode * rep_node = (*rep_chr).getLeftEnd();                        // make a point that starts at the centromere on the replicate chromosome

    while ( (*temp_node).getRightNode() != 0 )
    {
        CNode * next_node = new ( CNode::newAddyAssign() ) CNode( (*temp_node).getJunction() );  // copy the new chrNode to the replicate chromosome
        (*next_node).setLeftNode(rep_node);                              // link the new node to the preceeding replicate node
        (*rep_node).setRightNode(next_node);
        rep_node = (*rep_node).getRightNode();                            // advance the replicate Node, and the template Node
        temp_node = (*temp_node).getRightNode();
    }

    (*rep_node).setRightNode( (*rep_chr).getRightEnd() );                 // when at the telomere, link the last node replicated to the new telomere
    (*new_right).setLeftNode(rep_node);

    return rep_chr;
}

//  GENOTYPE A POSITION ON THE CHROMOSOME
int SexChromosome::positionAnc(double p)
{
    int anc = -1;

    // if it's an X or Z chromosome, check it like an autosome
    if ( m_type == X || m_type == Z )
    {
        CNode* right_node = (*m_left_end).getRightNode();

        while ( (*right_node).getJPosition() <= p && (*right_node).getRightNode() != 0 )
        {
            right_node = (*right_node).getRightNode();
        }

        CNode* left_node = (*right_node).getLeftNode();
        anc = (*left_node).getJAncestry();
    }
    else
    {
        if ( p <= s_par_boundary )
        {
            CNode* right_node = (*m_left_end).getRightNode();

            while ( (*right_node).getJPosition() <= p && (*right_node).getRightNode() != 0 )
            {
                right_node = (*right_node).getRightNode();
            }

            CNode* left_node = (*right_node).getLeftNode();
            anc = (*left_node).getJAncestry();
        }
        else
        {
            anc = (*m_right_end).getJAncestry() + s_num_anc;    // if it's beyond the par, it gets a whole different ancestry
        }
    }
    return anc;
}

// DISPLAY THE CHROMOSOME'S LIST OF JUNCTIONS IN THE TERMINAL WINDOW
void SexChromosome::displayChromosome()
{
    CNode *chr_node = m_left_end;
    Junction *jun;

    std::cout << "[" << m_type << ":" << m_number << "(" << m_length << "M)]_";

    while ( (*chr_node).getRightNode() != 0 )
    {
        jun = (*chr_node).getJunction();
        std::cout << "(" << (*jun).getChromosome() << "," << (*jun).getPosition() << "," << (*jun).getAncestry() << "," << (*jun).getNumOccur() <<  ")_" ;
        CNode *next_node = (*chr_node).getRightNode();
        chr_node = next_node;
    }

    jun = (*chr_node).getJunction();
    std::cout << "(" << (*jun).getChromosome() << "," << (*jun).getPosition() << "," << (*jun).getAncestry() << "," << (*jun).getNumOccur() << ")" << std::endl;
}

// PRINT THE CHROMOSOME'S LIST OF JUNCTIONS TO A .TXT FILE
void SexChromosome::textChromosome()
{
    CNode *chr_node = m_left_end;
    Junction *jun;

    g_check_file << "[" << m_type << ":" << m_number << "(" << m_length << "M)]_";

    while ( (*chr_node).getRightNode() != 0 )
    {
        jun = (*chr_node).getJunction();
        g_check_file << "(" << (*jun).getChromosome() << "," << (*jun).getPosition() << "," << (*jun).getAncestry() << "," << (*jun).getNumOccur() <<  ")_" ;
        CNode *next_node = (*chr_node).getRightNode();
        chr_node = next_node;
    }

    jun = (*chr_node).getJunction();
    g_check_file << "(" << (*jun).getChromosome() << "," << (*jun).getPosition() << "," << (*jun).getAncestry() << "," << (*jun).getNumOccur() << ")" << std::endl;
}

#endif // SEXCHROMOSOME_H
