#ifndef INTERACTION_GRAPHS_H
#define INTERACTION_GRAPHS_H

#include <map>
#include <algorithm>
#include <string>
#include <sstream>

#include "jungen_acc.cpp"
#include "gene.h"

/*
    INTGRAPH: interaction graph

    Interaction graphs are used to represent the effects of epistasis between paris of loci.
    The interaction graph holds information about the types of interactions and
    which genotypes have effects on a phenotype.

    Interactions are regarded as edges in a graph between potentially interacting alleles in
    a genotype. Depending on the model of epistasis, a proportion of a maximum effect is applied
    based on how many edges are present.
*/
class IntGraph
{
    public:
        // Interaction Graph Constructor prototype
        IntGraph( int p, Gene* locA, Gene* locB, int ancA, int ancB, double s, int d );

        // Setters

        // Getters
        double getMaxSel() { return m_max_sel; }
        Gene* getLocusA() { return m_locus_a; }
        Gene* getLocusB() { return m_locus_b; }
        int getAncA() { return m_anc_a; }
        int getAncB() { return m_anc_b; }
        double calc_Selection(int a1, int a2, int b1, int b2);
        int getPhenotype() { return m_phenotype; }
        int getEdgeFX() { return m_edge_fx; }

        // epistasis models
        double recEdges(int, double, int);  // returns the effect on the phenotype for a rec-rec interaction
        double domEdges(int, double, int);  // returns the effect ..                     dom-dom interaction
        double addEdges(int, double, int);  // returns the additive effets of the edges

        //static variables
        static int s_num_ancestries;


    protected:
	int m_phenotype;            // the phenotype the interaction impacts	
	Gene* m_locus_a;            // A pointer to the first gene in the pair
        Gene* m_locus_b;            // A pointer to the second gene in the pair	
	int m_anc_a;                // the ancestry at locus a that effects the interation
	int m_anc_b;                // the ancestry at the second locus ...        	
	double m_max_sel;           // The maximum allowed effect on the phenotype (usually fitness, hence max_sel)
        int m_edge_fx;              // the type of the effect of the interaction (edge)
        std::map<std::string, int> m_graph;   // a map between the genotypes and the number of edges they incur based on the model
        bool m_y_chromosome;        // whether or not the Y chromsome impacts the edge
};

// Iteraction Graph Constructor
IntGraph::IntGraph( int p, Gene* locA, Gene* locB, int a, int b, double s, int d ):
    m_phenotype(p), m_locus_a(locA), m_locus_b(locB), m_anc_a(a), m_anc_b(b), m_max_sel(s), m_edge_fx(d)
{

}


// calculate selection for key and genes given
double IntGraph::calc_Selection(int a1, int a2, int b1, int b2)
{
    int edges = 0;
    double value = 0;

    // increment the number of edges for each interacting pair of alleles in the genotype
    if( a1 == m_anc_a && b1 == m_anc_b )  { edges++; }
    if( a1 == m_anc_a && b2 == m_anc_b )  { edges++; }
    if( a2 == m_anc_a && b1 == m_anc_b )  { edges++; }
    if( a2 == m_anc_a && b2 == m_anc_b )  { edges++; }

    // apply the portion of selection according to the model of epistasis
    switch(m_edge_fx)
    {
        case(0):
            value = recEdges(edges, m_max_sel, 4);
            break;

        case(1):
            value = domEdges(edges, m_max_sel, 4);
            break;

        case(2):
            if ( edges > 1 && ( b1 == b2 || ( b1 >= s_num_ancestries || b2 >= s_num_ancestries ) ) )
            {
                value = m_max_sel;
            }
            break;

        case(3):
            value = addEdges(edges, m_max_sel, 4);
            break;

        default:
            value = 0;
            std::cout << "Edge Effect not properly specified!" << std::endl;
            break;
    }

    return value;
}

// Treat the effect of epistatic interactions as recessive
double IntGraph::recEdges(int e, double s, int f = 4)
{
    double v;
    if ( e == f ) { v = s;}
    else { v = 0; }
    return v;
}

// Treat the effect of epistatic interactons as dominant
double IntGraph::domEdges(int e, double s, int f = 4)
{
    double v;
    if ( e > 0 ) { v = s; }
    else { v = 0; }
    return v;
}

// Treat the effect of epistatic interactions as additive
double IntGraph::addEdges(int e, double s, int f = 4)
{
    double edges = e;
    double total = f;
    double v = s * ( edges / total );
    return v;
}

#endif // INTERACTION_GRAPHS_H
