
#ifndef CYTOPLASM_H
#define CYTOPLASM_H

#include "chromosome.h"

extern std::ofstream g_check_file;

class Cytoplasm : public Chromosome
{
    public:
        Cytoplasm(CNode* cen, CNode* tel, int* loc, int n, double l, ChrType t);
        Cytoplasm(Cytoplasm& a);
        virtual ~Cytoplasm();

        // getters
        double virtual getParB() { return 0.0; }
        double virtual getLength() { return 0.0; }
        ChrType virtual getType() {return m_type; }
        int virtual getNumber() { return m_number; }


        // display functions
        void virtual displayChromosome();
        void virtual textChromosome();

        // Operators
        //Cytoplasm& operator = ( Cytoplasm& c );
        virtual Chromosome* duplicateChr(int* l);
        virtual int positionAnc(double p);
        //virtual void recycle();


    private:

};

// CYTOPLASM CONSTRUCTOR
Cytoplasm::Cytoplasm(CNode * cen, CNode * tel, int* loc, int n, double l, ChrType t ):
    Chromosome(cen, tel, loc, n, 0, t)
{
    m_length = 0;
}

// CHROMOSOME COPY CONSTRUCTOR
Cytoplasm::Cytoplasm(Cytoplasm& a):
    Chromosome(a)
{

}

// CYTOPLASM DESTRUCTOR
Cytoplasm::~Cytoplasm()
{

}

// DUPLICATE A CYTOPLASM TO BE INHERITED
Chromosome* Cytoplasm::duplicateChr(int* l)
{
    // Get new end nodes
    CNode * new_left = new ( CNode::newAddyAssign() ) CNode( (*m_left_end).getJunction() );
    CNode * new_right = new ( CNode::newAddyAssign() ) CNode( (*m_right_end).getJunction() );

    // cytoplasmic genomes lack recombination and junctions
    // connect the two ends, and set the outside pointers to null
    (*new_left).setRightNode(new_right);
    (*new_left).setLeftNode(0);
    (*new_right).setLeftNode(new_left);
    (*new_right).setRightNode(0);
    Chromosome * rep_chr = new ( Chromosome::s_newAddyAssign() ) Cytoplasm(new_left, new_right, l, m_number, m_length, m_type);

    return rep_chr;
}

// GENOTYPE THE CYTOPLASM AT A PARTICULAR POSITION...
int Cytoplasm::positionAnc(double p)
{
    // no recombination, no genetic distance, just return the ancestry of the left end
    int a = (*m_left_end).getJAncestry();
    return a;
}

// DISPLAY THE CHROMOSOME'S LIST OF JUNCTIONS IN THE TERMINAL WINDOW
void Cytoplasm::displayChromosome()
{
    // only print the left and right end since those are the only two junctions
    CNode *chr_node = m_left_end;
    Junction *jun;

    std::cout << "[" << m_type << "_" << m_number << "]_";

    jun = (*m_left_end).getJunction();
    std::cout << "(" << (*jun).getChromosome() << "," << (*jun).getPosition() << "," << (*jun).getAncestry() << "," << (*jun).getNumOccur() <<  ")_" ;

    chr_node = m_right_end;
    jun = (*chr_node).getJunction();
    std::cout << "(" << (*jun).getChromosome() << "," << (*jun).getPosition() << "," << (*jun).getAncestry() << "," << (*jun).getNumOccur() << ")" << std::endl;
}

// PRINT THE CHROMOSOME'S LIST OF JUNCTIONS TO A .TXT FILE
void Cytoplasm::textChromosome()
{
    CNode *chr_node = m_left_end;
    Junction *jun;

    g_check_file << "[" << m_type << "_" << m_number << "]_";

    jun = (*m_left_end).getJunction();
    g_check_file << "(" << (*jun).getChromosome() << "," << (*jun).getPosition() << "," << (*jun).getAncestry() << "," << (*jun).getNumOccur() <<  ")_" ;

    chr_node = m_right_end;
    jun = (*chr_node).getJunction();
    g_check_file << "(" << (*jun).getChromosome() << "," << (*jun).getPosition() << "," << (*jun).getAncestry() << "," << (*jun).getNumOccur() << ")" << std::endl;
}

#endif //CYTOPLASM_H

