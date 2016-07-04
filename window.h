#ifndef WINDOW_H
#define WINDOW_H

#include "haphazenum.h"

// a window is used to return information about junction density on a segment of a chromosome
// it is not a part of the chromosome, but an object that is used to calculate junction density

struct Window
{
    Window (int c, ChrType t, double s, double e);
    int chromosome;
    ChrType cType;
    double start;
    double end;
    double numJunctions;
    double hyb_Index;
    int numOfChr;

    Window operator+(Window w);
    Window operator*(Window w);
    Window operator-(Window w);
    Window operator/(double n);
    int reset();
};

// WINDOW CONSTRUCTOR
Window::Window (int c, ChrType t, double s, double e) :
    chromosome(c), cType(t), start(s), end(e), numJunctions(0), hyb_Index(0), numOfChr(0)
{

}

// reset the windows counts to zero
int Window::reset()
{
    numJunctions = 0;
    hyb_Index = 0;
    numOfChr = 0;
}

// an addition operator for windows
Window Window::operator+(Window w)
{
    Window sum(w.chromosome, w.cType, w.start, w.end);
    if( chromosome == w.chromosome && start == w.start && end == w.end )
    {
        sum.numJunctions = numJunctions + w.numJunctions;
        sum.hyb_Index = hyb_Index + w.hyb_Index ;
        sum.numOfChr = numOfChr + w.numOfChr;
    }
    else
    {
        cout << "Attempt to add non-overlapping windows!" << endl;
    }
    return sum;
}

// an addition operator for windows
Window Window::operator*(Window w)
{
    Window product(w.chromosome, w.cType, w.start, w.end);
    if( chromosome == w.chromosome && start == w.start && end == w.end )
    {
        product.numJunctions = numJunctions * w.numJunctions;
        product.hyb_Index = hyb_Index * w.hyb_Index;
        product.numOfChr = numOfChr * w.numOfChr;
    }
    else
    {
        cout << "Attempt to multiply non-overlapping windows!" << endl;
    }
    return product;
}

// A subtraction operator for windows
Window Window::operator-(Window w)
{
    Window diff(w.chromosome, w.cType, w.start, w.end);
    if( chromosome == w.chromosome && start == w.start && end == w.end )
    {
        diff.numJunctions = numJunctions - w.numJunctions;
        diff.hyb_Index = hyb_Index - w.hyb_Index;
        diff.numOfChr = numOfChr - w.numOfChr;
    }
    else
    {
        cout << "Attempt to subtract non-overlapping windows!" << endl;
    }
    return diff;
}

// A division operator for windows
Window Window::operator/(double n)
{
    Window quotient(chromosome, cType, start, end);
    if( n != 0 )
    {
        quotient.numJunctions = numJunctions / n;
        quotient.hyb_Index = hyb_Index / n;
        quotient.numOfChr = numOfChr / n;
    }
    else
    {

        cout << "Attempt to divide window by zero!" << endl;
    }
    return quotient;
}

#endif //WINDOW_H
