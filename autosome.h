
#ifndef AUTOSOME_H
#define AUTOSOME_H

#include "chromosome.h"

extern std::ofstream g_check_file;

class Autosome : public Chromosome
{
    public:
        Autosome(CNode* cen, CNode* tel, int* loc, int n, double l, ChrType t);
        Autosome(Autosome& a);
        virtual ~Autosome();

        // getters
        int virtual getNumber() { return m_number; }
        double virtual getLength() { return m_length; }
        double virtual getParB() { return 1.0; }
        ChrType virtual getType() { return A; }

        // chromosome functions
        virtual int positionAnc(double p);


        // display functions
        void virtual displayChromosome();
        void virtual textChromosome();

        // Operators
        // Autosome& operator = ( Autosome& c );
        virtual Chromosome* duplicateChr(int* l);
        // virtual void recycle();


    private:

};

// AUTOSOME CONSTRUCTOR
Autosome::Autosome(CNode * cen, CNode * tel, int* loc, int n, double l, ChrType t = A):
    Chromosome(cen, tel, loc, n, l, t)
{

}

// CHROMOSOME DESTRUCTOR
Autosome::~Autosome()
{
  
}


// DUPLICATE AN AUTOSOME TO BE INHERITED
Chromosome* Autosome::duplicateChr(int *l)
{
    CNode * new_left = new ( CNode::newAddyAssign() ) CNode( (*m_left_end).getJunction() );          // make a new CNode that is a copy of the left end
    CNode * new_right = new ( CNode::newAddyAssign() ) CNode( (*m_right_end).getJunction()  );            // copy the right end
    Chromosome * rep_chr = new ( Chromosome::s_newAddyAssign() ) Autosome(new_left, new_right, l, m_number, m_length);       // make a new chromosome using the copied left and right end
    CNode * temp_node = (*m_left_end).getRightNode();                        // make pointer that starts on the CNode after the left end on the template chromosome
    CNode * rep_node = (*rep_chr).getLeftEnd();                        // make a point that starts at the left end on the replicate chromosome

    while ( (*temp_node).getRightNode() != 0 )
    {
        CNode * next_node = new (CNode::newAddyAssign() ) CNode( (*temp_node).getJunction() );  // copy the new chrNode to the replicate chromosome
        (*next_node).setLeftNode(rep_node);                              // link the new node to the preceeding replicate node
        (*rep_node).setRightNode(next_node);
        rep_node = (*rep_node).getRightNode();                            // advance the replicate Node, and the template Node
        temp_node = (*temp_node).getRightNode();
    }

    (*rep_node).setRightNode( (*rep_chr).getRightEnd() );                 // when at the right end, link the last node replicated to the new right end
    (*new_right).setLeftNode(rep_node);

    return rep_chr;
}

//  GENOTYPE A POSITION ON THE CHROMOSOME
int Autosome::positionAnc(double p)
{
    int anc = -1;

    CNode* right_node = (*m_left_end).getRightNode();

    while ( (*right_node).getJPosition() <= p && (*right_node).getRightNode() != 0 )
    {
        right_node = (*right_node).getRightNode();
    }

    CNode* left_node = (*right_node).getLeftNode();
    anc = (*left_node).getJAncestry();

    return anc;
}

// DISPLAY THE CHROMOSOME'S LIST OF JUNCTIONS IN THE TERMINAL WINDOW
void Autosome::displayChromosome()
{
    // start at the left
    CNode *chr_node = m_left_end;
    Junction *jun;

    // print a label for the chromosome
    std::cout << "[" << m_type << ":" << m_number << "(" << m_length << "M)]_";

    // move right printing info for each junction
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
void Autosome::textChromosome()
{
    // start at the left end
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

#endif //AUTOSOME_H
