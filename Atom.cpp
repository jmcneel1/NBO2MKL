#include "Atom.h"

using namespace std;

unsigned int Atom::GetMinS() const
{
  PeriodicTable ptable;
  return ptable.SCount(myElement);
}

unsigned int Atom::GetMinP() const
{
  PeriodicTable ptable;
  return ptable.PCount(myElement);
}

unsigned int Atom::GetMinD() const
{
  PeriodicTable ptable;
  return ptable.DCount(myElement);
}

unsigned int Atom::GetMinF() const
{
  PeriodicTable ptable;
  return ptable.FCount(myElement);
}

unsigned int Atom::GetMinCount() const
{
  PeriodicTable ptable;
  return ptable.ElecConfigCount(myElement);
}

string Atom::GetElement() const
{
  return myElement;
}

float Atom::GetX() const
{
  return myX;
}

float Atom::GetY() const
{
  return myY;
}

float Atom::GetZ() const
{
  return myZ;
}

float Atom::GetWeight() const
{
  return myWeight;
}

unsigned int Atom::GetAtomNumber() const
{
  return myAtomNumber;
}

string Atom::GetFullElement() const
{
  PeriodicTable ttable;
  return ttable.GetElement(myElement);
}

void Atom::ChangeX (float nx)
{
  myX = nx;
}

void Atom::ChangeY (float ny)
{
  myY = ny;
}

void Atom::ChangeZ (float nz)
{
  myZ = nz;
}

void Atom::ChangeElement (string newel)
{
  PeriodicTable ttable;
  if ( newel.length() > 2 )
  {
    myElement = ttable.GetElement(newel);
  }
  else
  {
    myElement = newel;
  }
  myAtomNumber = ttable.GetAtomNumber(myElement);
  myWeight = ttable.GetWeight(myAtomNumber);
}

void Atom::ChangeAtomNum (unsigned int newanum)
{
  PeriodicTable ttable;
  myAtomNumber = newanum;
  myElement = ttable.GetElement(newanum);
  myWeight = ttable.GetWeight(newanum);
}

Atom::Atom()
{
  myElement = "X";
  myAtomNumber = 0;
  myWeight = -1.000000;
  myX = 0;
  myY = 0;
  myZ = 0;
}

Atom::Atom(string input)
{
  PeriodicTable ttable;
  if ( input.length() > 2 )
  {
    myElement = ttable.GetElement(input);
  }
  else
  {
    myElement = input;
  }
  myAtomNumber = ttable.GetAtomNumber(myElement);
  myWeight = ttable.GetWeight(myAtomNumber);
  myX = 0;
  myY = 0;
  myZ = 0;
}

Atom::Atom(unsigned int anum)
{
  PeriodicTable ttable;
  myAtomNumber = anum;
  myElement=ttable.GetElement(anum);
  myWeight=ttable.GetWeight(anum);
  myX = myY = myZ = 0;
}

Atom::Atom(string el,float x, float y, float z)
{
  PeriodicTable ttable;
  if ( el.length() > 2 )
  {
    myElement = ttable.GetElement(el);
  }
  else
  {
    myElement = el;
  }
  myAtomNumber = ttable.GetAtomNumber(el);
  myWeight = ttable.GetWeight(myAtomNumber);
  myX = x;
  myY = y;
  myZ = z;
}

Atom::Atom(unsigned int anum, float x, float y, float z)
{
  PeriodicTable ttable;
  myAtomNumber = anum;
  myElement = ttable.GetElement(anum);
  myWeight = ttable.GetWeight(anum);
  myX = x;
  myY = y;
  myZ = z;
}

Atom::Atom(string el, float weight)
{
  PeriodicTable ttable;
  if ( el.length() > 2 )
  {
    myElement = ttable.GetElement(el);
  }
  else
  {
    myElement = el;
  }
  myAtomNumber = ttable.GetAtomNumber(myElement);
  myWeight = weight;
  myX = myY = myZ = 0;
}

Atom::Atom(unsigned int anum, float weight)
{
  PeriodicTable ttable;
  myAtomNumber = anum;
  myElement = ttable.GetElement(anum);
  myWeight = weight;
  myX = myY = myZ = 0;
}

Atom::Atom(string el, float x, float y, float z, float weight)
{
  PeriodicTable ttable;
  if ( el.length() > 2 )
  {
    myElement = ttable.GetElement(el);
  }
  else
  {
    myElement = el;
  }
  myAtomNumber = ttable.GetAtomNumber(myElement);
  myWeight = weight;
  myX = x;
  myY = y;
  myZ = z;
}

Atom::Atom(unsigned int anum, float x, float y, float z, float weight)
{
  PeriodicTable ttable;
  myAtomNumber = anum;
  myElement = ttable.GetElement(anum);
  myWeight = weight;
  myX = x;
  myY = y;
  myZ = z;
}

Atom& Atom::operator= (const Atom & rhs)
{
  myElement = rhs.myElement;
  myAtomNumber = rhs.myAtomNumber;
  myX = rhs.myX;
  myY = rhs.myY;
  myZ = rhs.myZ;
  myWeight = rhs.myWeight;
}
