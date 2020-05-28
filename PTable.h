#ifndef ORCANBO_PTABLE_H
#define ORCANBO_PTABLE_H

#include <tuple>
#include <string>
#include <vector>
#include <sstream>
#include <locale>

using namespace std;

class PeriodicTable
{
  public:
    float GetWeight(string el) const;
    float GetWeight(unsigned int atomnum) const;
    unsigned int GetAtomNumber(string el) const;
    string GetElement(string el) const;
    string GetElement(unsigned int atomnum) const;
    unsigned int ElecConfigCount(string el) const;
    unsigned int ElecConfigCount(unsigned int atomnum) const;
    unsigned int SCount(string el) const;
    unsigned int PCount(string el) const;
    unsigned int DCount(string el) const;
    unsigned int FCount(string el) const;
    unsigned int SCount(unsigned int atomnum) const;
    unsigned int PCount(unsigned int atomnum) const;
    unsigned int DCount(unsigned int atomnum) const;
    unsigned int FCount(unsigned int atomnum) const;
    PeriodicTable();
  private:
    vector<tuple<unsigned int,string,string,float> > myAtoms;
};

#endif
