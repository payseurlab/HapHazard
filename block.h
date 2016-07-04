#ifndef BLOCK_H
#define BLOCK_H

#include "haphazenum.h"

// blocks are used to return information about the block length distribution along the chromosome
struct Block
{
    int ancestry;
    int chromosome;
    void* haplotype;
    ChrType cType;
    double start;
    double end;
    void* startJunction;
    int startGen;
};

#endif //BLOCK_H
