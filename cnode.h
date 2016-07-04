
#ifndef CNODE003_H
#define CNODE003_H

#include <deque>
#include "junction.h"

// COMMENTS ON THE CNODE CLASS
//
// Background/Implementation notes:
// CNodes, chromosome nodes, are used to implement a chromosome as a doubly-linked list.
// They are modelled after the DLLs described in Drozdek 3rd. ed. with a couple of slight
// additions. Each CNode contains a pointer to the node in each direction on the chromosome: "left" and "right".
// The nodes that correspond to centromeres and telomeres have null
// pointers in the left and right orientation respectively. In addition, rather than containing their own
// data, a third pointer points to a junction in the common population "junction pool". By doing this, junctions
// can be used on multiple chromosomes, individuals, and sub-populations (demes) simultaneously. In addition,
// whenever a CNode is instantiated, it increments the number of occurences of the junction it points to, and
// decrements this variable when its destructor is called.

class CNode
{
  public:
    // constructor prototype
    CNode(Junction *j);
    // copy constructor prototype
    CNode(CNode& c);
    // destructor prototype
    ~CNode();

    // setters
    CNode* setRightNode ( CNode * d ) { m_right_node = d; return m_right_node; }
    CNode* setLeftNode ( CNode * p ) { m_left_node = p; return m_left_node; }
    Junction* setJunction (Junction * j);

    // getters
    Junction* getJunction() { return m_junction; }
    double getJPosition() { return (*m_junction).getPosition(); }
    double getJAncestry() {  return (*m_junction).getAncestry(); }
    CNode * getRightNode() { return m_right_node; }
    CNode * getLeftNode() { return m_left_node; }

    int locateJunction(int* l);


    static CNode* newAddyAssign();
    static void drainCNodeReservoir();

    static std::deque<CNode*> s_cnode_reservoir;
    static int s_num_cnode;
    void recycle();

  private:
    // variables
    Junction* m_junction;
    CNode* m_right_node;
    CNode* m_left_node;
    std::deque<int*>::iterator m_jun_location;
};



// CNODE CONSTRUCTOR
CNode::CNode(Junction *j):
    m_junction(j)
{
    //cout << this << " " << i << " " << location << endl;
    if (m_junction != 0 )
    {
        (*m_junction).incNumOccur();
    }

    m_left_node = 0;
    m_right_node = 0;
    s_num_cnode++;
}

// CNODE DESTRUCTOR -- makes sure to decrement the number of occurences of the junction
CNode::~CNode()
{
    (*m_junction).decNumOccur();

    s_num_cnode--;
}

// MAKE THE CNODE NULL -- POINT TO NOTHING -- BEFORE SENDING IT TO THE GARBAGE COLLECTOR
void CNode::recycle()
{
    if ( m_junction != 0 )
    {
        (*m_junction).decNumOccur();
        s_cnode_reservoir.push_back(this);
    }
    m_junction = 0;
    s_num_cnode--;
}

// GET NEXT PRE-ALLOCATED CNODE ADDRESS -- OR ALLOCATE NEW ONES AND RETURN THE LAST
CNode* CNode::newAddyAssign()
{
    CNode* new_addy;


    if ( ! s_cnode_reservoir.empty() )
    {
        new_addy = (*( s_cnode_reservoir.begin() ));
        s_cnode_reservoir.pop_front();
    }
    else
    {
        new_addy = new CNode(0);
    }

    return new_addy;
}

// CLEAN ALL CNODES IN RESERVED MEMORY
void CNode::drainCNodeReservoir()
{
    std::deque<CNode*>::iterator iter_c;
    for( iter_c = s_cnode_reservoir.begin() ; iter_c != s_cnode_reservoir.end() ; iter_c++ )
    {
        delete *iter_c;
    }
    s_cnode_reservoir.clear();
}

// JUNCTION SETTER FUNCTION -- if the junction is reset, we need to make sure the change in occurence is accounted for
Junction* CNode::setJunction (Junction * j)
{
    (*m_junction).decNumOccur();
    m_junction = j;
    if ( j != 0 )
    {
        (*m_junction).incNumOccur();
    }
    return m_junction;
}

int CNode::locateJunction(int* l)
{
    return (*m_junction).locate(l);
}
#endif // CNODE003_H
