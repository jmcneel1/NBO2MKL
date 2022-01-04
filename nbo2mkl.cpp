#include "BasisSet.h"
#include "Atom.h"
#include <sstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <cstdlib>

using namespace std;

// This writes the header for the MKL file
// Pretty self-explanatory
void WriteHeader ( ofstream & outFile )
{
  outFile << "$MKL" << endl;
  outFile << "#" << endl;
  outFile << "# MKL format from an NBO Plot File" << endl;
  outFile << "#" << endl;
}

/*
This procedure takes arguments to the function
to define the charge and multiplicity that
are written to the output file (the .mkl file)
*/

void WriteChargeMult ( ofstream & outFile, char* charge, char* mult )
{
  string scharge(charge);
  string smult(mult);
  outFile << "$CHAR_MULT" << endl;
  outFile << "  "+scharge+" "+smult << endl;
  outFile << "$END" << endl << endl;
}

/*
This reads the .47 file for the coordinates
It assumes that they are in Angstroems
It then write them into the output file (.mkl)
and saves the geometry to a vector
*/

void ReadWriteGeom ( char* base, ofstream & outFile,
                     vector<Atom> & geom )
{
  string ifile(base);
  ifile = ifile+".47";
  ifstream inFile(ifile.c_str());
  string line;
  outFile << "$COORD" << endl;
  // The GoTo mautil function, then read 2 more lines
  GoTo(inFile,"$COORD");
  getline(inFile,line); getline(inFile,line);
  while ( line.find("$END") == string::npos )
  {
    unsigned int atnum, wt;
    float x, y, z;
    stringstream ss;
    ss << line;
    ss >> atnum >> wt >> x >> y >> z;
    Atom tatom(atnum,x,y,z);
    geom.push_back(tatom);
    outFile << setw(4) << atnum;
    outFile << setw(12) << fixed << setprecision(6) << x;
    outFile << setw(11) << fixed << setprecision(6) << y;
    outFile << setw(11) << fixed << setprecision(6) << z << endl;
    getline(inFile,line);
  }
  outFile << "$END" << endl << endl;
  inFile.close();
}

/*
The following section includes ReadWriteBasis
And all the corresponding procedures used in the
function
*/

/*
Here we read the number of basis functions from the NBAS entry
in the $GENNBO keylist
*/

unsigned int ReadNBAS ( ifstream & inFile )
{
  string temp("");
  while ( temp.find("NBAS") == string::npos )
  {
    inFile >> temp;
  }
  size_t pos = temp.find("=");
  pos++;
  return atoi(temp.substr(pos).c_str());
}

/*
Here we read the centers for the archive file... there should
be nbasis centers to read
*/

void ReadCenters( ifstream & inFile, vector<unsigned int> & centers,
                  unsigned int nbas)
{
  unsigned int center;
  GoTo(inFile,"$BASIS");
  string temp1, temp2;
  // should be temp1 == "CENTER" temp2 == "="
  inFile >> temp1 >> temp2;
  for ( unsigned int i = 0; i < nbas; i++ )
  {
    inFile >> center;
    centers.at(i) = center;
  }
}

/*
Now we read the basis set labels and also return the
max angular momentum
*/

void ReadLabels ( ifstream & inFile, vector<unsigned int> & labels,
                  unsigned int & maxlabel, unsigned int nbas )
{
  unsigned int label;
  maxlabel = 0; // Just in case...
  string temp1, temp2;
  // should be temp1 == "LABEL" temp2 == "="
  inFile >> temp1 >> temp2;
  for ( unsigned int i = 0; i < nbas; i++ )
  {
    inFile >> label;
    if ( label/100 > maxlabel ) maxlabel = label/100;
    labels.at(i) = label;
  }
}

/*
This function assigns the correct shell to each BF
Makes things easier later on...
We make the following assumptions which is based on the
way ORCA prints out it's MOs
 - s shells begin/end with 1
 - p shells begin with 103 and end with 102
 - d shells begin with 255 and end with 251
 - high shells go 1...N in sequential order
*/

void AssignShells ( vector<unsigned int> & shell,
                    const vector<unsigned int> & labels )
{
  unsigned int currentshell(-1);
  for ( unsigned int i = 0; i < labels.size(); i++ )
  {
    if ( labels.at(i) == 1 || labels.at(i) == 103 ||
         labels.at(i) == 255 || labels.at(i) == 351 ||
         labels.at(i) == 451 || labels.at(i) == 551 ) currentshell++;
    shell.at(i) = currentshell;
  }
}

/*
Here we read the NSHELL variable from the archive
*/

void ReadNShell ( ifstream & inFile, unsigned int & nshell )
{
  GoTo(inFile,"$CONTRACT");
  string t1, t2;
  // This should be t1 == NSHELL, t2 == =
  inFile >> t1 >> t2 >> nshell;
}

/*
Here we read the NEXP variable from the archive
*/

void ReadNExp ( ifstream & inFile, unsigned int & nexp)
{
  string t1, t2;
  //This should be t1 == NEXP, t2 == =
  inFile >> t1 >> t2 >> nexp;
}

/*
Here we read the number of primitives of the shells
*/

void ReadNPrim ( ifstream & inFile, vector<unsigned int> & nprim )
{
  unsigned int temp;
  string t1, t2;
  // First we skip the NCOMP entries
  inFile >> t1 >> t2;
  for ( unsigned int i = 0; i < nprim.size(); i++ ) inFile >> temp;
  // Now we get there...
  inFile >> t1 >> t2;
  for ( unsigned int i = 0; i < nprim.size(); i++ )
  {
    inFile >> temp;
    nprim.at(i) = temp;
  }
}

/*
Here we read the ptr of the shell to the relevant
primitive set
*/

void ReadPTR ( ifstream & inFile, vector<unsigned int> & ptr )
{
  unsigned int temp;
  string t1, t2;
  inFile >> t1 >> t2;
  for ( unsigned int i = 0; i < ptr.size(); i++ )
  {
    inFile >> temp;
    ptr.at(i) = temp;
  }
}

/*
Now we actually read the exponents
*/

void ReadExp ( ifstream & inFile,
               vector<double> & exps, unsigned int nexp)
{
  string t1, t2;
  double temp;
  // Should be t1 == EXP, t2 == =
  inFile >> t1 >> t2;
  for ( unsigned int i = 0; i < nexp; i++ )
  {
    inFile >> temp;
    exps.at(i) = temp;
  }
}

/*
This procedure will read the contraction coefficients for
one of the angular momentum sets ... only stores the coefficient
if it's nonzero. This is checked by comparing absolute values of
the existing value with the read value
*/
void ReadShellCoeff ( ifstream & inFile,
                      vector<double> & coeffs, unsigned int nexp)
{
  double temp;
  string t1, t2;
  // Should be t1 == CL, t2 == =
  inFile >> t1 >> t2;
  for ( unsigned int i = 0; i < nexp; i++ )
  {
    inFile >> temp;
    if ( fabs(coeffs.at(i)) < fabs(temp) ) coeffs.at(i) = temp;
  }
}

/*
Here we actually organize all of our information and create
the basis set object
*/

void SetBasisSetObject( BasisSet & basis, const vector<double> & exps,
                        const vector<double> & coeffs,
                        const vector<unsigned int> & centers,
                        const vector<unsigned int> & labels,
                        const vector<unsigned int> & shell,
                        const vector<unsigned int> & ptr,
                        const vector<unsigned int> & nprim )
{
  for ( unsigned int i = 0; i < labels.size(); i++ )
  {
    vector<double> texp(nprim.at(shell.at(i)));
    vector<double> tcoeff(nprim.at(shell.at(i)));
    for ( unsigned int j = 0; j < texp.size(); j++ )
    {
      texp.at(j) = exps.at(ptr.at(shell.at(i))-1+j);
      tcoeff.at(j) = coeffs.at(ptr.at(shell.at(i))-1+j);
    }
    basis.AddBF(labels.at(i),centers.at(i),texp,tcoeff);
  }
}

/*
This reads the basis set information from the .47
file and stores it in the BasisSet object
--- This assumes spherical basis sets (ORCA)
*/

void ReadBasis ( ifstream & inFile,
                 const vector<Atom> & geom,
                 BasisSet & basis )
{
  // We assume here that NBAS is present in the archive
  // This is required by GENNBO unless REUSE=T, so
  // we assume that REUSE=F here ...
  unsigned int nbas = ReadNBAS(inFile);
  vector<unsigned int> centers(nbas); // length of nbasis
  ReadCenters(inFile,centers,nbas);
  unsigned int maxL; // Stores the highest L
  vector<unsigned int> labels(nbas); // length of nbasis
  ReadLabels(inFile,labels,maxL,nbas);
  vector<unsigned int> shell(nbas);
  AssignShells(shell,labels);
  unsigned int nshell, nexp;
  ReadNShell(inFile,nshell);
  ReadNExp(inFile,nexp);
  vector<unsigned int> nprim(nshell);
  ReadNPrim(inFile,nprim);
  vector<unsigned int> ptr(nshell);
  ReadPTR(inFile,ptr);
  vector<double> exps(nexp,0);
  ReadExp(inFile,exps,nexp);
  vector<double> coeffs(nexp,0);
  for ( unsigned int i = 0; i <= maxL; i++ )
  {
    ReadShellCoeff(inFile,coeffs,nexp);
  }
  SetBasisSetObject(basis,exps,coeffs,centers,labels,shell,ptr,nprim);
}

/*
Now we write the basis to the .mkl file
We only write based on the label == 1, 103, 255, 351, 451, 551
*/

void WriteBasis ( ofstream & outFile, const vector<Atom> & geom,
                  const BasisSet & basis )
{
  unsigned int crt_center = basis.GetCenter(0);
  outFile << "$BASIS" << endl;
  for ( unsigned int i = 0; i < basis.GetNBasis(); i++ )
  {
    if ( basis.GetCenter(i) != crt_center )
    {
      crt_center = basis.GetCenter(i);
      outFile << "$$\n";
    }
    if ( basis.GetLabel(i) == 1 )
    {
      outFile << " 1 S 1.0" << endl;
      for ( unsigned int j = 0; j < basis.GetNumPrim(basis.GetShell(i)); j++ )
      {
        outFile << setw(19) << fixed << setprecision(9) << basis.GetExponent(basis.GetPtr(basis.GetShell(i))+j);
        outFile << setw(17) << fixed << setprecision(9) << basis.GetNormContractionS(basis.GetPtr(basis.GetShell(i))+j) << endl;
      }
    }
    else if ( basis.GetLabel(i) == 103 )
    {
      outFile << " 3 P 1.0" << endl;
      for ( unsigned int j = 0; j < basis.GetNumPrim(basis.GetShell(i)); j++ )
      {
        outFile << setw(19) << fixed << setprecision(9) << basis.GetExponent(basis.GetPtr(basis.GetShell(i))+j);
        outFile << setw(17) << fixed << setprecision(9) << basis.GetNormContractionP(basis.GetPtr(basis.GetShell(i))+j) << endl;
      }
    }
    else if ( basis.GetLabel(i) == 255 )
    {
      outFile << " 5 D 1.0" << endl;
      for ( unsigned int j = 0; j < basis.GetNumPrim(basis.GetShell(i)); j++ )
      {
        outFile << setw(19) << fixed << setprecision(9) << basis.GetExponent(basis.GetPtr(basis.GetShell(i))+j);
        outFile << setw(17) << fixed << setprecision(9) << basis.GetNormContractionD(basis.GetPtr(basis.GetShell(i))+j) << endl;
      }
    }
    else if ( basis.GetLabel(i) == 351 )
    {
      outFile << " 7 F 1.0" << endl;
      for ( unsigned int j = 0; j < basis.GetNumPrim(basis.GetShell(i)); j++ )
      {
        outFile << setw(19) << fixed << setprecision(9) << basis.GetExponent(basis.GetPtr(basis.GetShell(i))+j);
        outFile << setw(17) << fixed << setprecision(9) << basis.GetNormContractionF(basis.GetPtr(basis.GetShell(i))+j) << endl;
      }
    }
    else if ( basis.GetLabel(i) == 451 )
    {
      outFile << " 9 G 1.0" << endl;
      for ( unsigned int j = 0; j < basis.GetNumPrim(basis.GetShell(i)); j++ )
      {
        outFile << setw(19) << fixed << setprecision(9) << basis.GetExponent(basis.GetPtr(basis.GetShell(i))+j);
        outFile << setw(17) << fixed << setprecision(9) << basis.GetNormContractionG(basis.GetPtr(basis.GetShell(i))+j) << endl;
      }
    }
    else if ( basis.GetLabel(i) == 551 )
    {
      outFile << "11 H 1.0" << endl;
      for ( unsigned int j = 0; j < basis.GetNumPrim(basis.GetShell(i)); j++ )
      {
        outFile << setw(19) << fixed << setprecision(9) << basis.GetExponent(basis.GetPtr(basis.GetShell(i))+j);
        outFile << setw(17) << fixed << setprecision(9) << basis.GetNormContractionH(basis.GetPtr(basis.GetShell(i))+j) << endl;
      }
    }
  }
  outFile << "$END" << endl << endl;
}

/*
  The main ReadWriteBasis function simply reads the .47 file
  into a BasisSet object and then calls a function
  to write it to the output file
*/

void ReadWriteBasis ( char* base, ofstream & outFile,
                      const vector<Atom> & geom,
                      BasisSet & basis )
{
  string ifile(base);
  ifile = ifile+".47";
  ifstream inFile(ifile.c_str());
  ReadBasis(inFile,geom,basis);
  WriteBasis(outFile,geom,basis);
  inFile.close();
}

/*
This routine actually reads the orbitals from the relevant file.
We then back-convert a few of the coefficients for compatibility
with ORCA
Namely -- we rearrange the order of p-orbitals and we
alter the sign of a couple of coefficients.
Remember ... ORCA wants the MKL file to contain x, y, z
But they print out z, x, y, and that is how they are stored in NBO
Thus we have z + 2, x -1
*/

void ReadOrbs ( ifstream & inFile,
                vector< vector<double> >  & orbs,
                const BasisSet & basis)
{
  // First we skip the first three lines
  string line;
  for ( unsigned int i = 0; i < 3; i++ ) getline(inFile,line);
  // Now define the temp double variable
  double temp;
  // Now we read the coefficients and alter relevant ones
  // based on the associated labels
  for ( unsigned int i = 0; i < basis.GetNBasis(); i++ )
  {
    for ( unsigned int j = 0; j < basis.GetNBasis(); j++ )
    {
      inFile >> temp;
      if ( basis.GetLabel(j) == 357 ) temp*=-1;
      if ( basis.GetLabel(j) == 356 ) temp*=-1;
      if ( basis.GetLabel(j) == 456 ) temp*=-1;
      if ( basis.GetLabel(j) == 457 ) temp*=-1;
      if ( basis.GetLabel(j) == 458 ) temp*=-1;
      if ( basis.GetLabel(j) == 459 ) temp*=-1;
      if ( basis.GetLabel(j) == 101 ) orbs.at(i).at(j-1) = temp;
      else if ( basis.GetLabel(j) == 102 ) orbs.at(i).at(j-1) = temp;
      else if ( basis.GetLabel(j) == 103 ) orbs.at(i).at(j+2) = temp;
      else orbs.at(i).at(j) = temp;
    }
  }
}

/*
And finally we write the orbitals to the mkl file
*/

void WriteOrbs ( ofstream & outFile, vector< vector<double > > & orbs,
                 unsigned int nbasis )
{
  // First the header line
  outFile << "$COEFF_ALPHA" << endl;
  // The first loop goes over the 'full' rows ( 5 columns )
  for ( unsigned int i = 0; i < nbasis/5; i++ )
  {
    // this might be sloppy ... but it works ... the formatting that is
    outFile << " ";
    for ( unsigned int j = 0; j < 5; j++ ) outFile << "a1g  ";
    outFile << endl;
    // Here we set each orbital as having 0 hartree energy
    for ( unsigned int j = 0; j < 5; j++ ) outFile << setw(13) << fixed << setprecision(7) << 0.0 << " ";
    outFile << " " << endl;
    for ( unsigned int j = 0; j < nbasis; j++ )
    {
      outFile << setw(12) << fixed << setprecision(7) << orbs.at(i*5).at(j);
      outFile << setw(13) << fixed << setprecision(7) << orbs.at(i*5+1).at(j);
      outFile << setw(13) << fixed << setprecision(7) << orbs.at(i*5+2).at(j);
      outFile << setw(13) << fixed << setprecision(7) << orbs.at(i*5+3).at(j);
      outFile << setw(13) << fixed << setprecision(7) << orbs.at(i*5+4).at(j);
      outFile << " " << endl;
    }
  }
  if ( nbasis%5 != 0 )
  {
    outFile << " ";
    for ( unsigned int i = 0; i < nbasis%5; i++ ) outFile << "a1g  ";
    outFile << endl;
    for ( unsigned int i = 0; i < nbasis%5; i++ ) outFile << setw(13) << fixed << setprecision(7) << 0.0 << " ";
    outFile << " " << endl;
    for ( unsigned int i = 0; i < nbasis; i++ )
    {
      outFile << setw(12) << fixed << setprecision(7) << orbs.at((nbasis/5)*5).at(i);
      for ( unsigned int j = 1; j < nbasis%5; j++ )
      {
        outFile << setw(13) << fixed << setprecision(7) << orbs.at((nbasis/5)*5+j).at(i);
      }
      outFile << endl;
    }
  }
  outFile << " $END" << endl;
}

/*
This is the routine to read the orbital coefficients
and write them to the mkl file
*/

void ReadWriteOrbitals ( char* base, ofstream & outFile,
                        const vector<Atom> & geom,
                        const BasisSet & basis,
                        char* ext )
{
  string ifile(base), ifileext(ext), line;
  ifile = ifile+"."+ifileext;
  ifstream inFile(ifile);
  vector< vector<double> > orbs(basis.GetNBasis(), vector<double> (basis.GetNBasis(),0.0));
  ReadOrbs(inFile,orbs,basis);
  WriteOrbs(outFile,orbs,basis.GetNBasis());
  inFile.close();
}

/*
This is the main routine.
This program is used to take plot files
generated by GENNBO (or an NBO interface
in common ESS packages) and create an MKL
file that can be used by ORCA's
orca_2mkl routine to generate a corresponding
gbw file that can be used as a guess for
for a quantum calculation.

In particular, this program was written to
utilize NBOs and NLMOs as guess orbitals
for CASSCF calculations.

Building: This was successfully built with gcc 7.4.0
  - Please feel free to try other versions!

Calling:
nbo2mkl BASE EXT charge multiplicity
  - Here BASE is the base of the relevant files.
  - It is required that there be the following present:
    1. BASE.EXT
    2. BASE.47
  - EXT is the extension of the relevant orbital file.
    - By default, GENNBO uses 37 for NBO and 39
      for NLMOs
  - the BASE.47 file is the archive file that was either
    used with GENNBO to generate the plot files or almost
    all ESS packages have an option to save the archive
    file
  - charge is the charge of the system
  - multiplicity is the multiplicity

Currently, the MKL file ONLY  contains the coordinates,
charge, multiplicity, basis set, and MO coefficients.

Also, it doesn't include symmetry labels or orbital
energies

Also, it doesn't contain ECP support.

The inclusion of the above along with occupancies is planned
for future development.
*/



int main ( int argc, char* argv[] )
{
  if ( argc != 5 )
  {
    cout << "Incorrect argument #!!!\nExiting..." << endl;
    exit(1);
  }
  BasisSet basis;
  vector<Atom> geom;
  string ofile(argv[1]);
  ofile = ofile + ".mkl";
  ofstream outFile(ofile.c_str());
  WriteHeader(outFile);
  WriteChargeMult(outFile,argv[2],argv[3]);
  ReadWriteGeom(argv[1],outFile,geom);
  ReadWriteBasis(argv[1],outFile,geom,basis);
  ReadWriteOrbitals(argv[1],outFile,geom,basis,argv[4]);
  outFile.close();
}
