
#ifndef JUNCTION003_H
#define JUNCTION003_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <deque>

// declare global variables
extern int g_generation;
extern std::ofstream g_check_file;

// COMMENTS ON THE JUNCTION CLASS
//
// Background notes:
// The junction is the most fundamental and important datatype in this program.
// A junction is a switchpoint between two different ancestries that occurs along the length
// of a chromosome. They are the result of recombination between chromosomal tracts of different
// ancestry, and are inherited as point mutations. As such, they provide a direct measure of the
// amount of hybridization between two diverged lineages of organisms.
//
// Implementation notes:
// The junction is defined by the following simple, numeric parameters:
// 1) Position in the genome
//     a) The chromsome it lies on and b) its genetic (or physical) position on that chromosome
// 2) Ancestry - the ancestry that follows the "right" of the junction on the chromsome
//
// The above three parameters are required to instantiate a junction object.
//
// In addition, we keep track of other parameters
// generation -- the generation the junction was born in
// m_num_occur -- the number of times the junction occurs in the (meta-)population
//
// Being a fundamental object, it is simple, and few functions beyond getters and setters are needed

// JUNCTION CLASS DEFINITION
class Junction
{
    public:
        // junction constructor prototype
        Junction(int c, double p, int a);
        ~Junction();

        // setters
        std::deque<Junction*>::iterator setPoolPosition( std::deque<Junction*>::iterator j ) { m_pool_position = j; return m_pool_position; }
        int setChromosome(int c) { m_chromosome = c; return m_chromosome; }
        double setPosition(double p) { m_position = p; return m_position; }
        int setAncestry(int a) { m_ancestry = a; return m_ancestry; }
        int setGen(int g) { m_junction_generation = g; return m_junction_generation; }

        // getters
        int getChromosome() { return m_chromosome; }
        double getPosition() { return m_position; }
        int getAncestry() { return m_ancestry; }
        int getGen() { return m_junction_generation; }
        int getNumOccur() { return m_num_occur; }
        int* getLocationCounts() { return m_location_counts; }
        std::deque<Junction*>::iterator getPoolPosition() { return m_pool_position; }

        // functions
        int incNumOccur() { m_num_occur++; return m_num_occur; }  // increase the count of this junction
        int decNumOccur(); // decrease the count of this unciton
        int incLocCount(int* i) { m_location_counts[*i]++; return m_location_counts[*i]; } // increase the count of this junction in a specific location
        int decLocCount(int* i) { m_location_counts[*i]--; return m_location_counts[*i]; } // decrease the count of this junciton in a specific location
        void displayJunction();     // print the junction's info on a terminal
        void textJunction();        // print the junction's info in a file
        int deallocCN(int* l);           // signal from a CNode that the junction has been deallocated
        void resetLocations();      // reset the location counts to zero for this junction
        int locate(int* l) { m_location_counts[*l]++; return m_location_counts[*l]; }       // set the location/location counts for the junction

        //STATIC junction functions and variables
        static int s_num_locations;             // the number of demes
        void resetLocationCounts();    // reset location counts for ALL junctions

        //Junction Pool functions
        // The junction pool holds all junctions that are in use by the simulation
        // The pool allows access to junction information without the need to access them
        // through individuals and chromosomes
        static void s_cleanJunctionPool();     // remove unused junctions from the pool
        static void s_drainJunctionPool();     // remove all junctions from the pool
        static void s_textJunctionPool();      // text all the junctions to the check_file

        // Junction Pool variables
        static std::deque<Junction*> s_junction_pool;        // the junction pool
        static int s_num_junctions;                     // the total number of junctions in the simulation

        // Junction recycling functions
        static Junction* s_newAddyAssign();             // returns the next available junction slot or creates a new one
        static void s_drainJunctionReservoir();                      // removes all junctions from the reservoir

        // Junction recycling variables
        static std::deque<Junction*> s_junction_reservoir;   // the reserve of junction memory slots
        static int s_num_recycled;                      // the number of junctions that were recycled

    private:
        int m_chromosome;           // the chromosome on which the junction resides
        double m_position;          // the genetic position of the junction
        int m_ancestry;             // the ancestry of the tract following the junction
        int m_junction_generation;  // the generation in which the junction was born
        int m_location;
        int m_num_occur;             // the number of times the junction is used in a simulation
        int* m_location_counts;      // an array containing the number of times the junction occurs in each deme
        std::deque<Junction*>::iterator m_pool_position;     // an iterator that points to a junction in the
                                                        // junction pool vector
};

// JUNCTION CONSTRUCTOR
Junction::Junction(int c, double p, int a):
    m_chromosome(c),
    m_position(p),
    m_ancestry(a),
    m_junction_generation(g_generation)
{
    m_num_occur = 0;   // initialize the number of occurrences to zero

    m_location_counts = new int[s_num_locations];      // create an array of junction counts for each deme
    resetLocations();                        // and set the location counts equal to zero
    s_junction_pool.push_back(this);   // add the new junction to the junction pool
    setPoolPosition(s_junction_pool.end() - 1); // get the iterator for this junction in the pool
    s_num_junctions++; // and increase the total count of junctions
}

// JUNCTION DESTRUCTOR
Junction::~Junction()
{
    s_num_junctions--; // decrease the total number of junctions
}

// DEALLOCATE A JUNCTION FROM A CNODE -- CALLED IN CNODE::DEALLOCATE()
int Junction::deallocCN(int* l)
{
    // decrease overall count and by location
    m_location_counts[*l]--;
    decNumOccur();
    return m_num_occur;
}

// RESET THE LOCATION COUNTS TO ZERO FOR ALL JUNCTIONS
void Junction::resetLocationCounts()
{
    std::deque<Junction*>::iterator iterJ;
    for( iterJ = s_junction_pool.begin() ; iterJ != s_junction_pool.end() ; iterJ++ )
    {
        (**iterJ).resetLocations();
    }
}

// RESET THE LOCATION COUNTS TO ZERO FOR A SINGLE JUNCTION
void Junction::resetLocations()
{
    for( int i = 0 ; i < (s_num_locations) ; i++ )
    {
        m_location_counts[i] = 0;
    }
}

// DECREASE THE TOTAL NUMBER OF OCCURENCES FOR THE JUNCTION AND RECYCLE IT IF NECESSARY
int Junction::decNumOccur()
{
    m_num_occur--;
    if(m_num_occur == 0 )
    {
        s_junction_reservoir.push_back(this);
        Junction::s_num_recycled++;
        s_num_junctions--;
    }

    return m_num_occur;
}

// GET NEXT PRE-ALLOCATED JUNCTION ADDRESS -- OR ALLOCATE NEW ONES AND RETURN THE LAST
Junction* Junction::s_newAddyAssign()
{
    Junction* new_addy = 0; // initialize a null junciton address pointer

    if ( s_junction_reservoir.size() != 0 )    // if there are empty junctions reserved
    {
        new_addy = (*( s_junction_reservoir.begin() )); // get the new address from the back of the vector
        s_junction_reservoir.pop_front(); // then remove it from the reservoir
    }
    else
    {

        new_addy = new Junction(-1,-1,-1); // otherwise, make a new junction with nonsense parameters
    }

    return new_addy;
}

// DISPLAY A JUNCTION IN THE TERMINAL
void Junction::displayJunction()
{
    std::cout << "(" << m_chromosome << ", " << m_position << ", " << m_ancestry << ", " << m_num_occur << ")" << std::endl;
}

// DISPLAY A JUNCTION IN THE TERMINAL
void Junction::textJunction()
{
    g_check_file << "\t" << m_chromosome << "\t" << m_position << "\t" << m_ancestry << "\t" << m_junction_generation << "\t" << m_num_occur;
    for( int i = 0 ; i < s_num_locations ; i++ )
    {
        g_check_file << "\t" << m_location_counts[i];
    }
    g_check_file << std::endl;
}

// JUNCTION POOL FUNCTIONS
// REMOVE EXTINCT JUNCTIONS FROM THE JUNCTION POOL
void Junction::s_cleanJunctionPool()
{
    // iterate over the pool and remove all those whose m_num_occur == 0
    int pool_size = s_junction_pool.size();
    std::deque<Junction*>::iterator iter_j;
    for ( iter_j = s_junction_pool.begin() ; iter_j != s_junction_pool.end() ; iter_j++ )
    {
        while ( (*(*iter_j)).getNumOccur() == 0 )
        {
            s_junction_reservoir.push_back(*iter_j);
            iter_j = s_junction_pool.erase(iter_j);
        }
    }
    std::cout << "JUNCTION POOL CLEANED " << std::endl;
    std::cout << "Start: " << pool_size << "\tNew Pool Size: " << s_junction_pool.size() << std::endl;
}

// REMOVE ALL JUNCTIONS FROM MEMORY
void Junction::s_drainJunctionPool()
{
    std::deque<Junction*>::iterator iter_j;
    int j = s_junction_pool.size();
    std::cout << s_junction_pool.size() << std::endl;
    int i = 0;
    for (iter_j = s_junction_pool.begin(); iter_j != s_junction_pool.end() ; iter_j++ )
    {
        std::cout << i << " of " << j << " cleared. " << std::endl;
        delete (*iter_j);
        i++;
    }
    s_junction_pool.clear();
}

// REMOVE ALL JUNCTIONS FROM MEMORY
void Junction::s_drainJunctionReservoir()
{
    std::deque<Junction*>::iterator iter_j;
    for (iter_j = s_junction_reservoir.begin(); iter_j != s_junction_reservoir.end() ; iter_j++ )
    {
        delete (*iter_j);
    }
    s_junction_reservoir.clear();
}

// PRINT ALL THE JUNCTIONS TO A TEXT FILE
void Junction::s_textJunctionPool()
{
    g_check_file << "Junction Pool " << " " << g_generation << std::endl;
    std::deque<Junction*>::iterator iter_j;
    for (iter_j = s_junction_pool.begin(); iter_j < s_junction_pool.end() ; iter_j++ )
    {
        g_check_file << (*iter_j) << "\t" ;
        (*(*iter_j)).textJunction();
    }
}

#endif // JUNCTION003_H
