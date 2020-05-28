#ifndef ORCANBO_ATOM_H
#define ORCANBO_ATOM_H

#include "PTable.h"

using namespace std;

class Atom
{
  public:
    // Return the element symbol from private variable myElement
    string GetElement() const;
    // Return X-coordinate
    float GetX() const;
    // Return Y-coordinate
    float GetY() const;
    // Return Z-coordinate
    float GetZ() const;
    // Return weight of nucleus
    float GetWeight() const;
    // Return Atomic Number
    unsigned int GetAtomNumber() const;
    // Return Full Name of Atom ... All UpperCase
    string GetFullElement() const;
    // Change X-coordinate
    void ChangeX (float nx);
    // Change Y-coordinate
    void ChangeY (float ny);
    // Change Z-Coordinate
    void ChangeZ (float nz);
    // Change Element .. this string can be full or symbol.
    // Elemnt is change, atom number is changed, weight is set to default
    void ChangeElement (string newel);
    // Change nuclar weight
    void ChangeWeight(float weight);
    // Change Atom Num ... also changes Element and Weight to default
    void ChangeAtomNum(unsigned int newanum);
    // This returns the minimum number of basis functions needed to describe
    // the ground state of the atom
    // H - 1, He -1, Li - 2, Be - 2, B - 5, etc...
    unsigned int GetMinCount () const;
    unsigned int GetMinS () const;
    unsigned int GetMinP () const;
    unsigned int GetMinD () const;
    unsigned int GetMinF () const;
    Atom();
    Atom(string input);
    Atom(unsigned int anum);
    Atom(string el,float x, float y, float z);
    Atom(unsigned int anum, float x, float y, float z);
    Atom(string el, float weight);
    Atom(unsigned int anum, float weight);
    Atom(string el, float x, float y, float z, float weight);
    Atom(unsigned int anum, float x, float y, float z, float weight);
    Atom& operator= (const Atom & rhs);
  private:
    string myElement;
    float myX,myY,myZ;
    unsigned int myAtomNumber;
    float myWeight;
};

#endif
