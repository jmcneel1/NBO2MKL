#ifndef BASISSET_NBO2MKL_H
#define BASISSET_NBO2MKL_H

#include <vector>
#include <iostream>
#include "mautil.h"
#include <math.h>
#include <stdlib.h>

using namespace std;

class BasisSet
{
  public:
    // This will calculate the normalization coefficient for
    // a single primitive with a given exponent (Gaussian exp(-a*r^2))
    // and a given angular momentum determined from label
    double Norm ( double exponent, unsigned int label ) const;
    // This is another version that calculated the norm for
    // 2 different exponential factors on the same center.
    double Norm ( double exp1, double exp2, unsigned int shell ) const;
    // This will return the normalization for a contracted basis function
    double NormContracted ( unsigned int shell,
                            const vector<double> & coeffs,
                            const vector<double> & exps ) const;
    // Returns the number of shells in the object
    unsigned int GetNShell() const { return myNshell; }
    // Returns the shell associated with basis function index
    unsigned int GetShell( unsigned int index ) const;
    // Returns the number of exponents (exp(-a*r^2)) stored
    // Equal to the total number of primitives (sum(shell*nprim(shell)))
    unsigned int GetNExp() const { return myNexp; }
    // Returns the number of basis functions.
    // Equals the length of the center and label vectors
    unsigned int GetNBasis() const { return myNBasis; }
    // returns the number of primitives for a
    // given shell.
    // return 20,000 if not good
    unsigned int GetNumPrim(unsigned int index) const;
    // If index is good, returns the number of components for a shell.
    // Is equal to 2*L+1
    unsigned int GetNumComponent(unsigned int index) const;
    // If index is good, Atomic center associated with basis function.
    // Counting starts at 0
    unsigned int GetCenter(unsigned int index) const;
    // If index is good, returns the normalized primitive coefficient.
    // Only actually normalized if not contracted...
    double GetNormContractionS ( unsigned int index ) const;
    // If index is good, Returns the Contraction given by ORCA...
    // In general, these aren't normalized, but we check anyway
    // in case the user input a normalized basis set...
    double GetContractionS ( unsigned int index ) const;
    double GetNormContractionP ( unsigned int index ) const;
    double GetContractionP ( unsigned int index ) const;
    double GetNormContractionD ( unsigned int index ) const;
    double GetContractionD ( unsigned int index ) const;
    double GetNormContractionF ( unsigned int index ) const;
    double GetContractionF ( unsigned int index ) const;
    double GetNormContractionG ( unsigned int index ) const;
    double GetContractionG ( unsigned int index ) const;
    double GetNormContractionH ( unsigned int index ) const;
    double GetContractionH ( unsigned int index ) const;
    // If index is good, returns alpha for a primitive
    double GetExponent(unsigned int index) const;
    // If index is good, returns the label for a basis function
    unsigned int GetLabel(unsigned int index) const;
    // Gets the ptr to the primitive that begins for the shell
    // first checks that index is valid
    unsigned int GetPtr(unsigned int index) const;
    // The following add shells, NOT basis functions
    void AddSShell(unsigned int center,
              const vector<double> & exps,
              const vector<double> & coeffs);
    void AddPShell(unsigned int center,
              const vector<double> & exps,
              const vector<double> & coeffs);
    void AddDShell(unsigned int center,
              const vector<double> & exps,
              const vector<double> & coeffs);
    void AddFShell(unsigned int center,
              const vector<double> & exps,
              const vector<double> & coeffs);
    void AddGShell(unsigned int center,
              const vector<double> & exps,
              const vector<double> & coeffs);
    void AddHShell(unsigned int center,
              const vector<double> & exps,
              const vector<double> & coeffs);
    // Adds a basis function
    // We assume that the shell components are ordered in the same way that
    // ORCA orders them:
    // P: 103, 101, 102
    // D: 255, 253, 252, 254, 251
    // F: 351, 352, 353, 354, 355, 356, 357
    // G: 451, 452, 453, 454, 455, 456, 457, 458, 459
    // H: 551, 552, 553, 554, 555, 556, 557, 558, 559, 560, 561
    void AddBF (unsigned int label, unsigned int center,
                const vector<double> & exps, const vector<double> & coeffs );
    // Constructors
    BasisSet ();
    BasisSet ( const BasisSet & obas );
    BasisSet & operator=(const BasisSet & obas );
  private:
    unsigned int myNBasis;
    unsigned int myNshell;
    unsigned int myNexp;
    vector<unsigned int> myNprim; // 1 for each shell
    vector<unsigned int> myNcomp; // 1 for each shell
    vector<unsigned int> myCenter; // 1 for each BF
    vector<unsigned int> myPtr; // 1 for each shell
    vector<double> myContractS; // 1 for each shell
    vector<double> myNContractS; // 1 for each shell
    vector<double> myContractP; // 1 for each shell
    vector<double> myNContractP; // 1 for each shell
    vector<double> myContractD; // 1 for each shell
    vector<double> myNContractD; // 1 for each shell
    vector<double> myContractF; // 1 for each shell
    vector<double> myNContractF; // 1 for each shell
    vector<double> myContractG; // 1 for each shell
    vector<double> myNContractG; // 1 for each shell
    vector<double> myContractH; // 1 for each shell
    vector<double> myNContractH; // 1 for each shell
    vector<double> myExponent; // 1 for each shell
    vector<unsigned int> myLabel; // 1 for each basis function
    vector<unsigned int> myShells; // Counting starts at 0
};

#endif
