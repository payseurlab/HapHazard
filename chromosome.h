
#ifndef CHROMOSOME003_H
#define CHROMOSOME003_H

#include <vector>
#include <algorithm>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "junction.h"
#include "cnode.h"

#include "block.h"
#include "haphazenum.h"

// COMMENTS ON THE CHROMOSOME CLASS
//
// Background notes:
// The chromosome in these simuations provides a means and structure for storing and manipulation
// information about ancestry throughout the genome. As such, they are thought of as ordered lists
// of junctions. Using junctions, we can ignore intevening genotypic information, thus making the
// chromsomes much smaller than they might be in other types of forward simulations. They also provide
// a scaffold upon which recombination is to occur, as well as being the main structure that is
// inherited by an individual's offspring
//
// This program will also feature special types of chromosomes: the sex chromosomes as well as cytoplasmic
// inheritance mechanism. Sex chromosomes will contain loci that indicate the sex of individuals, and will
// also not recombine in heterogametic sexes. Cytoplasmic factors will be inherited in single copy, either
// maternally or paternally depending on the biology of the organism being studied
//
// Implementation notes:
// Each chromosome is implemented (roughly speaking) as a circular, doubly-linked list. The chromosome object,
// itself, is not this list but merely the starting point for operating on the list. Three parameters are required
// to instantiate an chromosome object. They are:
//      1) centromere - a pointer to a CNode which indicates the 0 M position on a chromosome(or chromosome arm--
//                note: the program was developed with mouse genomes in mind, of which all chromosomes are
//                      telocentric. This will be adapted in later versions to accomodate other chromosomal
//                      architectures -- metacentric, etc.
//      2) telomere - a pointer to CNode which indicates the distal end of the chromosome
//      3) type - a variable of type ChrType ( an enum defined in jgenum003.h) indicating the type of chromosome
//              i.e. autosomal vs. sex chromosome or cytoplasmic. The type may also imply the length.
//      4) length - the genetic length of the chromosome in Morgans, also the position of the telomere
//
// Of the above, 1, 2, and 3 are required to instantiate a chromosome. 4 is implied by three.
//
// Chromosomes also provide functions to calculate basic summary statistics, such as number of junctions, and
// the average tract length of ancestry tracts on the chromosome. They can also make use of markers used for
// interrogating ancestry at specific positions in the genome. In addition, chromosomes provide functions for
// displaying the information on them.
//
//

extern std::ofstream g_check_file;


class Chromosome
{
    public:
        // Chromosome constructor prototype
        Chromosome(CNode * left, CNode * right, int* loc, int n, double l, ChrType t);

        // Chromosome destructor prototype
        virtual ~Chromosome();

        // Setters
        CNode* setLeftEnd(CNode* c) { m_left_end = c; return m_left_end; }
        CNode* setRightEnd(CNode* t) { m_right_end = t; return m_right_end; }

        // Getters
        CNode* getLeftEnd() { return m_left_end; }
        CNode* getRightEnd() { return m_right_end; }
        virtual double getLength() { return -1;}
        virtual ChrType getType() { return A; }
        virtual int getNumber() { return -1; }
        virtual double getParB() { return -1; }
        int* getLocation() { return m_location;}

        // Biological functions
        virtual Chromosome* duplicateChr(int* l) { return 0; }
        virtual int positionAnc(double p) { return -1; }  // return the ancestry at a position p (any position)

        // Summary Statistic Calcuators
        int calcNumJunctions();
        double calcATL();
        std::vector<Block> getBlocks(int c);
        void locateJunctions();

        // Data management and display
        virtual void displayChromosome() {};
        virtual void textChromosome() {};
        std::vector<int> junByWindow(double w);
        void deleteNodes();
        void recycle();

        //STATIC
        static void s_drainChromosomeReservoir();

        static std::deque<Chromosome*> s_chromosome_reservoir;
        static Chromosome* s_newAddyAssign();
        static double s_par_boundary;
        static gsl_rng* s_rand_num;
        static double s_window_size;
        static int s_num_anc;

        // Operators
        //Chromosome& operator = (Chromosome& c);


    protected:
        CNode* m_left_end;      // pointer to the left end node
        CNode* m_right_end;     // pointer to the right end node
        int m_number;           // the genomic index
        double m_length;        // the length of the chromosome
	ChrType m_type;         // the type of the chromosome (X, Y, A etc.)	
	int* m_location;        // the deme THIS chromosome resides in
};

// CHROMOSOME CONSTRUCTOR
Chromosome::Chromosome(CNode * cen, CNode * tel, int* loc, int n, double l, ChrType t):
    m_left_end(cen), m_right_end(tel), m_number(n), m_length(l), m_type(t), m_location(loc)
{

}

// CHROMOSOME DESTRUCTOR
Chromosome::~Chromosome()
{
    // Iterate over the nodes from left to right and recycle them
    CNode* chr_node = m_left_end;

    while ( (*chr_node).getRightNode() != 0 ) // the rightmost node has a null right pointer
    {
        CNode* dump_node = chr_node;
        chr_node = (*dump_node).getRightNode();
        (*dump_node).recycle();
    }

    (*chr_node).recycle();
}

// RECYCLE THE CHROMOSOME
void Chromosome::recycle()
{
    // iterate over the nodes from left to right, and recycle each node
    CNode* chrNode = m_left_end;
    while ( (*chrNode).getRightNode() != 0 ) // the rightmost node has a null right pointer
    {
        CNode *dumpNode = chrNode;
        chrNode = (*dumpNode).getRightNode();
        (*dumpNode).recycle();
    }
    (*chrNode).recycle();

    // when done, set the ends to null
    (*this).setLeftEnd(0);
    (*this).setRightEnd(0);

    // and add this one to the reservoir
    s_chromosome_reservoir.push_back(this);
}

// REMOVE ALL THE NODES ON THE CHROMOSOME FROM MEMORY
void Chromosome::deleteNodes()
{
    CNode *chr_node = m_left_end;

    // iterate left to right, call the destructor to delete nodes
    while ( (*chr_node).getRightNode() != 0 ) // the rightmost node has a null right pointer
    {
        CNode *del_node = chr_node;
        chr_node = (*del_node).getRightNode();
        delete del_node;
    }

    delete chr_node;
}

// CLEAN ALL CHROMOSOMES FROM RESERVED MEMORY
void Chromosome::s_drainChromosomeReservoir()
{
    // iterator over the reservoir deque, and call delete for each chromosome pointed to
    std::deque<Chromosome*>::iterator iter_c;
    for( iter_c = s_chromosome_reservoir.begin() ; iter_c != s_chromosome_reservoir.end() ; iter_c++ )
    {
        delete *iter_c;
    }

    // after the actual chromosomes are gone, clear the reservoir
    s_chromosome_reservoir.clear();
}

// GET A CHROMOSOME BLOCK FROM THE CHROMOSOME RESERVOIR
Chromosome* Chromosome::s_newAddyAssign()
{
    Chromosome* new_addy;

    // check the reservoir
    if ( s_chromosome_reservoir.size() != 0 ) // if it's not empty
    {
        new_addy = (*( s_chromosome_reservoir.begin() )); // get the first chromosome block
        s_chromosome_reservoir.pop_front(); // and remove it from the reservoir
    }
    else    // otherwise
    {
        new_addy = new Chromosome(0,0,0,-1,-1,A);   // make a new chromosome
    }

    return new_addy;
}

// UPDATE THE LOCATION COUNTS FOR THE JUNCTIONS ON THIS CHROMOSOME
void Chromosome::locateJunctions()
{
    CNode* node = m_left_end;
    while( (*node).getRightNode() != 0 )
    {
        (*node).locateJunction(m_location);
        node = (*node).getRightNode();
    }

    (*node).locateJunction(m_location);
}

// COUNT THE NUMBER OF JUNCTIONS ON THE CHROMOSOME
int Chromosome::calcNumJunctions()
{
    int n = 0; // start the count at zero
    CNode* chr_node = m_left_end; // start counting at the left end

    while ( (*chr_node).getRightNode() != 0 ) // keep moving right until you reach the rightmost node
    {
        if ( (*chr_node).getJunction() != 0 ) // if the node has a junction
        {
            n++;                                // count it
        }
        CNode* next_node = (*chr_node).getRightNode();
        chr_node = next_node;
    }

    return n;
}

// COUNT THE NUMBER OF JUNCTION PER WINDOW ON THE CHROMOSOME
std::vector<int> Chromosome::junByWindow(double w)
{
    std::vector<int> windows; // initialize a vector to store windows
    CNode * node = m_left_end;
    int i = 1;

    // iterate over windows by counting, and adding their lengths until the total is >= the
    // length of the chromosome
    for ( i = 1 ; i * w < m_length ; i++ )
    {
        int n = 0; // initialize the count for this window to zero
        do
        {
            node = (*node).getRightNode();
            if ( (*node).getJPosition() < i * w ) // check to see if the junction is in the window
            {
                n++; // if so, count it
            }
            // keep going until the position of the node is outside the window
        } while ( (*node).getJPosition() <= i * w && (*node).getRightNode() != 0 );

        windows.push_back(n);
    }

    return windows;
}


// GET THE BLOCKS ON THE CHROMOSOME
std::vector<Block> Chromosome::getBlocks(int c)
{
    // each block resides between two junctions
    // define the first for the first two junctions
    // and assign data to the block
    std::vector<Block> blocks;

    CNode* start = m_left_end;
    CNode* end = (*m_left_end).getRightNode();


    Block newBlock;
    newBlock.ancestry = (*start).getJAncestry();
    newBlock.chromosome = m_number;
    newBlock.haplotype = this;
    newBlock.cType = m_type;
    newBlock.start = (*start).getJPosition();
    newBlock.startJunction = (*start).getJunction();
    newBlock.startGen = (*( (*start).getJunction() )).getGen();
    newBlock.end = (*end).getJPosition();
    blocks.push_back(newBlock);
    start = end;

    // iterate over the rest of the chromosome assigning block data
     while ( (*end).getRightNode() != 0 )
    {
        end = (*end).getRightNode();
        newBlock.ancestry = (*start).getJAncestry();
        newBlock.chromosome = m_number;
        newBlock.haplotype = this;
        newBlock.cType = m_type;
        newBlock.startJunction = (*start).getJunction();
        newBlock.startGen = (*( (*start).getJunction() )).getGen();
        newBlock.start = (*start).getJPosition();
        newBlock.end = (*end).getJPosition();

        blocks.push_back(newBlock);
        start = end;
    }

    return blocks;
}

// CALCULATE THE AVERAGE ANCESTRY TRACT LENGTH ON THE CHROMOSOME
double Chromosome::calcATL()
{
    int n = calcNumJunctions();
    double atl = m_length / n;
    return atl;
}


#endif // CHROMOSOME003_H
