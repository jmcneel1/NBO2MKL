#include <iostream>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>

using namespace std;

void WriteChargeMult ( ofstream & outFile, char* charge, char* mult )
{
  string scharge(charge);
  string smult(mult);
  outFile << "$CHAR_MULT" << endl;
  outFile << "  "+scharge+" "+smult << endl;
  outFile << "$END" << endl << endl;
}

void WriteHeader ( ofstream & outFile )
{
  outFile << "$MKL" << endl;
  outFile << "#" << endl;
  outFile << "# MKL format from an NBO Plot File" << endl;
  outFile << "#" << endl;
}

void ReadWriteGeom ( char* base, ofstream & outFile, unsigned int & numatom )
{
  numatom = 0;
  string ifile(base);
  ifile = ifile+".31";
  ifstream inFile(ifile.c_str());
  string line;
  outFile << "$COORD" << endl;
  for ( unsigned int i = 0; i < 6; i++ ) { getline(inFile,line); }
  while ( line.find("---") == string::npos )
  {
    numatom++;
    unsigned int atnum;
    double x, y, z;
    stringstream ss;
    ss << line;
    ss >> atnum >> x >> y >> z;
    outFile << setw(4) << atnum;
    outFile << setw(12) << fixed << setprecision(6) << x;
    outFile << setw(11) << fixed << setprecision(6) << y;
    outFile << setw(11) << fixed << setprecision(6) << z << endl;
    getline(inFile,line);
  }
  outFile << "$END" << endl << endl;
  inFile.close();
}

void ReadWriteCharge ( char* base, ofstream & outFile, unsigned int natom )
{
  string ifile(base);
  ifile = ifile+".out";
  ifstream inFile(ifile.c_str());
  string line;
  outFile << "$CHARGES" << endl;
  getline(inFile,line);
  while ( line.find("Atom No    Charge        Core      Valence    Rydberg      Total") == string::npos )
  {
    getline(inFile,line);
  }
  getline(inFile,line);
  for ( unsigned int i = 0; i < natom; i++ )
  {
    string at; unsigned int index; double charge;
    getline(inFile,line);
    stringstream ss;
    ss << line;
    ss >> at >> index >> charge;
    outFile << setw(11) << fixed << setprecision(6) << charge << endl;
  }
  outFile << "$END" << endl << endl;
  inFile.close();
}

void ReadWriteBasis ( char* base, ofstream & outFile, unsigned int natom,
                      vector<unsigned int> & labels )
{
  string ifile(base);
  vector<unsigned int> aindex, nprim, shell;
  vector<double> exps, coeffs;
  unsigned int maxshell(0), numprim(0);
  ifile = ifile+".31";
  ifstream inFile(ifile.c_str());
  string line;
  getline(inFile,line);
  while (line.find("---") == string::npos) { getline(inFile,line); }
  getline(inFile,line);
  while	(line.find("---") == string::npos) { getline(inFile,line); }
  getline(inFile,line);
  while (line.find("---") == string::npos) { getline(inFile,line); }
  getline(inFile,line);
  while ( line.find("---") == string::npos )
  {
    unsigned int center, ncom, ptr, prim, label;
    stringstream ss;
    ss << line;
    ss >> center >> ncom >> ptr >> prim;
    aindex.push_back(center);
    nprim.push_back(prim);
    numprim+=prim;
    getline(inFile,line);
    ss.clear(); ss.str("");
    ss << line;
    ss >> label;
    shell.push_back(label/100);
    if ( maxshell < label/100 ) maxshell = label/100;
    if ( label == 1 ) labels.push_back(1);
    else if ( label == 103 )
    {
      labels.push_back(103); labels.push_back(101); labels.push_back(102);
    }
    else if ( label == 255 )
    {
      labels.push_back(255); labels.push_back(253); labels.push_back(252);
      labels.push_back(254); labels.push_back(251);
    }
    else if ( label == 351 )
    {
      labels.push_back(351); labels.push_back(352); labels.push_back(353);
      labels.push_back(354); labels.push_back(355); labels.push_back(356);
      labels.push_back(357);
    }
    else if ( label == 451 )
    {
      labels.push_back(451); labels.push_back(452); labels.push_back(453);
      labels.push_back(454); labels.push_back(455); labels.push_back(456);
      labels.push_back(457); labels.push_back(458); labels.push_back(459);
    }
    else if ( label == 551 )
    {
      labels.push_back(551); labels.push_back(552); labels.push_back(553);
      labels.push_back(554); labels.push_back(555); labels.push_back(556);
      labels.push_back(557); labels.push_back(558); labels.push_back(559);
      labels.push_back(560); labels.push_back(561);
    }
    getline(inFile,line);
  }
  for ( unsigned int i = 0; i < numprim/4; i++ )
  {
    getline(inFile,line);
    stringstream ss;
    double t1, t2, t3, t4;
    ss << line;
    ss >> t1 >> t2 >> t3 >> t4;
    exps.push_back(t1); exps.push_back(t2);
    exps.push_back(t3); exps.push_back(t4);
  }
  if ( numprim%4 != 0 )
  {
    getline(inFile,line);
    stringstream ss;
    ss << line;
    for ( unsigned int i = 0; i < numprim%4; i++ )
    {
      double t1;
      ss >> t1;
      exps.push_back(t1);
    }
  }
  getline(inFile,line);
  for ( unsigned int j = 0; j < numprim/4; j++ )
  {
    getline(inFile,line);
    stringstream ss;
    double t1, t2, t3, t4;
    ss << line;
    ss >> t1 >> t2 >> t3 >> t4;
    coeffs.push_back(t1); coeffs.push_back(t2);
    coeffs.push_back(t3); coeffs.push_back(t4);
  }
  if ( numprim%4 != 0 )
  {
    getline(inFile,line);
    stringstream ss;
    ss << line;
    for ( unsigned int i = 0; i < numprim%4; i++ )
    {
      double t1;
      ss >> t1;
      coeffs.push_back(t1);
    }
  }
  getline(inFile,line);
  for ( unsigned int i = 1; i <= maxshell; i++ )
  {
    for ( unsigned int j = 0; j < numprim/4; j++ )
    {
      getline(inFile,line);
      stringstream ss;
      double t1, t2, t3, t4;
      ss << line;
      ss >> t1 >> t2 >> t3 >> t4;
      if ( abs(coeffs.at(j*4)) < abs(t1) ) coeffs.at(j*4) = t1;
      if ( abs(coeffs.at(j*4+1)) < abs(t2) ) coeffs.at(j*4+1) = t2;
      if ( abs(coeffs.at(j*4+2)) < abs(t3) ) coeffs.at(j*4+2) = t3;
      if ( abs(coeffs.at(j*4+3)) < abs(t4) ) coeffs.at(j*4+3) = t4;
    }
    if ( numprim%4 != 0 )
    {
      getline(inFile,line);
      stringstream ss;
      ss << line;
      for ( unsigned int j = 0; j < numprim%4; j++ )
      {
        double t1;
        ss >> t1;
        if ( abs(coeffs.at((numprim/4)*4+j)) < abs(t1) ) coeffs.at((numprim/4)*4+j) = t1;
      }
    }
    getline(inFile,line);
  }
  inFile.close();
  outFile << "$BASIS" << endl;
  unsigned int ccenter = aindex.at(0);
  unsigned int count = 0, primcount(0);
  while ( count < aindex.size() )
  {
    if ( ccenter == aindex.at(count) )
    {
      if ( shell.at(count) == 0 ) outFile << " 1 S 1.0" << endl;
      else if ( shell.at(count) == 1 ) outFile << " 3 P 1.0" << endl;
      else if ( shell.at(count) == 2 ) outFile << " 5 D 1.0" << endl;
      else if ( shell.at(count) == 3 ) outFile << " 7 F 1.0" << endl;
      else if ( shell.at(count) == 4 ) outFile << " 9 G 1.0" << endl;
      else if ( shell.at(count) == 5 ) outFile << "11 H 1.0" << endl;
      for ( unsigned int i = 0; i < nprim.at(count); i++ )
      {
        outFile << setw(19) << fixed << setprecision(9) << exps.at(primcount);
        outFile << setw(17) << fixed << setprecision(9) << coeffs.at(primcount) << endl;
        primcount++;
      }
    }
    else
    {
      outFile << "$$" << endl;
      if ( shell.at(count) == 0 ) outFile << " 1 S 1.0" << endl;
      else if ( shell.at(count) == 1 ) outFile << " 3 P 1.0" << endl;
      else if ( shell.at(count) == 2 ) outFile << " 5 D 1.0" << endl;
      else if ( shell.at(count) == 3 ) outFile << " 7 F 1.0" << endl;
      else if ( shell.at(count) == 4 ) outFile << " 9 G 1.0" << endl;
      else if ( shell.at(count) == 5 ) outFile << "11 H 1.0" <<	endl;
      for ( unsigned int i = 0; i < nprim.at(count); i++ )
      {
        outFile << setw(19) << fixed << setprecision(9) << exps.at(primcount);
        outFile << setw(17) << fixed << setprecision(9) << coeffs.at(primcount) << endl;
        primcount++;
      }
      ccenter++;
    }
    count++;
  }
  outFile << "$END" << endl << endl;
}

void ReadWriteOrbitals ( char* base, ofstream & outFile, unsigned int natom,
                         const vector<unsigned int> & labels )
{
  string t1, t2, t3;
  string ifile(base);
  ifile = ifile+".47";
  ifstream inFile(ifile.c_str());
  string line;
  getline(inFile,line);
  stringstream ss; ss << line;
  ss >> t1 >> t2 >> t3;
  unsigned int nbasis =  atoi(t3.substr(5).c_str());
  inFile.close();
  vector< vector<double> > nbos(nbasis, vector<double> (nbasis,0.0));
  string ifile2(base);
  unsigned int orb(0);
  unsigned int row(0);
  ifile2 = ifile2+".37";
  inFile.open(ifile2.c_str());
  for ( unsigned int i = 0; i < 3; i++ ) getline(inFile,line);
  for ( unsigned int i = 0; i < nbasis; i++ )
  {
    for ( unsigned int j = 0; j < nbasis; j++ )
    {

      double temp;
      inFile >> temp;
      if ( labels.at(j) == 357 ) temp*=-1;
      if ( labels.at(j) == 356 ) temp*=-1;
      if ( labels.at(j) == 456 ) temp*=-1;
      if ( labels.at(j) == 457 ) temp*=-1;
      if ( labels.at(j) == 458 ) temp*=-1;
      if ( labels.at(j) == 459 ) temp*=-1;
      if ( labels.at(j) == 101 ) nbos.at(i).at(j-1) = temp;
      else if ( labels.at(j) == 102 ) nbos.at(i).at(j-1) = temp;
      else if ( labels.at(j) == 103 ) nbos.at(i).at(j+2) = temp;
      else nbos.at(i).at(j) = temp;
    }
  }
  inFile.close();
  outFile << "$COEFF_ALPHA" << endl;
  for ( unsigned int i = 0; i < nbasis/5; i++ )
  {
    outFile << " ";
    for ( unsigned int j = 0; j < 5; j++ ) outFile << "a1g  ";
    outFile << endl;
    for ( unsigned int j = 0; j < 5; j++ ) outFile << setw(13) << fixed << setprecision(7) << 0.0 << " ";
    outFile << " " << endl;
    for ( unsigned int j = 0; j < nbasis; j++ )
    {
      outFile << setw(12) << fixed << setprecision(7) << nbos.at(i*5).at(j);
      outFile << setw(13) << fixed << setprecision(7) << nbos.at(i*5+1).at(j);
      outFile << setw(13) << fixed << setprecision(7) << nbos.at(i*5+2).at(j);
      outFile << setw(13) << fixed << setprecision(7) << nbos.at(i*5+3).at(j);
      outFile << setw(13) << fixed << setprecision(7) << nbos.at(i*5+4).at(j);
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
  }
  outFile << " $END" << endl;
}

int main ( int argc, char* argv[] )
{
  unsigned int natom;
  vector<unsigned int> labels;
  string ofile(argv[1]);
  ofile = ofile + ".mkl";
  ofstream outFile(ofile.c_str());
  WriteHeader(outFile);
  WriteChargeMult(outFile,argv[2],argv[3]);
  ReadWriteGeom(argv[1],outFile,natom);
  ReadWriteCharge(argv[1],outFile,natom);
  ReadWriteBasis(argv[1],outFile,natom,labels);
  ReadWriteOrbitals(argv[1],outFile,natom,labels);
  return 0;
}
