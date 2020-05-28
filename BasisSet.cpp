#include "BasisSet.h"

double BasisSet::Norm ( double exponent, unsigned int label ) const
{
  unsigned int l = label/100;
  double a = pow(2*exponent/3.14159265359,0.75);
  double b = pow(pow(2,l)/DFac(2*l-1),0.5);
  double c = pow(2*exponent,l/2.0);
  return (a*b*c);
}

double BasisSet::Norm ( double exp1, double exp2, unsigned int shell) const
{
  double a = pow((exp1+exp2)/3.14159265359,0.75);
  double b = pow(pow(2,shell)/DFac(2*shell-1),0.5);
  double c = pow((exp1+exp2),shell/2.0);
  return(a*b*c);
}

double BasisSet::NormContracted ( unsigned int shell,
                        const vector<double> & coeffs,
                        const vector<double> & exps ) const
{
  double total = 0;
  for ( unsigned int i = 0; i < coeffs.size(); i++ )
  {
    for ( unsigned int j = 0; j < coeffs.size(); j++ )
    {
      total+=coeffs.at(i)*coeffs.at(j)*pow(1/Norm(exps.at(i),exps.at(j),shell),2);
    }
  }
  return total;
}

unsigned int BasisSet::GetNumPrim(unsigned int index) const
{
  if ( index < myNshell )
  {
    return myNprim.at(index);
  }
  else
  {
    cout << "Invalid Index For GenNumPrim!" << endl << "Exiting..." << endl;
    exit(1);
  }
}

unsigned int BasisSet::GetNumComponent(unsigned int index) const
{
  if ( index < myNshell )
  {
    return myNcomp.at(index);
  }
  else
  {
    cout << "Invalid Index For GenNumComponent!" << endl << "Exiting..." << endl;
    exit(1);
  }
}

unsigned int BasisSet::GetCenter(unsigned int index) const
{
  if ( index < myNBasis )
  {
    return myCenter.at(index);
  }
  else
  {
    cout << "Invalid index for GetCenter!" << endl << "Exiting..." << endl;
    exit(1);
  }
}

double BasisSet::GetNormContractionS ( unsigned int index ) const
{
  if ( index < myNexp )
  {
    return myNContractS.at(index);
  }
  else
  {
    cout << "Invalid index for GetNormContractionS!" << endl << "Exiting...\n";
    exit(1);
  }
}

double BasisSet::GetContractionS ( unsigned int index ) const
{
  if ( index < myNexp )
  {
    return myContractS.at(index);
  }
  else
  {
    cout << "Invalid index for GetContracttionS!" << endl << "Exiting...\n";
    exit(1);
  }
}

double BasisSet::GetNormContractionP ( unsigned int index ) const
{
  if ( index < myNexp )
  {
    return myNContractP.at(index);
  }
  else
  {
    cout << "Invalid index for GetNormContractionP!" << endl << "Exiting...\n";
    exit(1);
  }
}

double BasisSet::GetContractionP ( unsigned int index ) const
{
  if ( index < myNexp )
  {
    return myContractP.at(index);
  }
  else
  {
    cout << "Invalid index for GetContracttionP!" << endl << "Exiting...\n";
    exit(1);
  }
}

double BasisSet::GetNormContractionD ( unsigned int index ) const
{
  if ( index < myNexp )
  {
    return myNContractD.at(index);
  }
  else
  {
    cout << "Invalid index for GetNormContractionD!" << endl << "Exiting...\n";
    exit(1);
  }
}

double BasisSet::GetContractionD ( unsigned int index ) const
{
  if ( index < myNexp )
  {
    return myContractD.at(index);
  }
  else
  {
    cout << "Invalid index for GetContracttionD!" << endl << "Exiting...\n";
    exit(1);
  }
}

double BasisSet::GetNormContractionF ( unsigned int index ) const
{
  if ( index < myNexp )
  {
    return myNContractF.at(index);
  }
  else
  {
    cout << "Invalid index for GetNormContractionF!" << endl << "Exiting...\n";
    exit(1);
  }
}

double BasisSet::GetContractionF ( unsigned int index ) const
{
  if ( index < myNexp )
  {
    return myContractF.at(index);
  }
  else
  {
    cout << "Invalid index for GetContracttionF!" << endl << "Exiting...\n";
    exit(1);
  }
}

double BasisSet::GetNormContractionG ( unsigned int index ) const
{
  if ( index < myNexp )
  {
    return myNContractG.at(index);
  }
  else
  {
    cout << "Invalid index for GetNormContractionG!" << endl << "Exiting...\n";
    exit(1);
  }
}

double BasisSet::GetContractionG ( unsigned int index ) const
{
  if ( index < myNexp )
  {
    return myContractG.at(index);
  }
  else
  {
    cout << "Invalid index for GetContracttionG!" << endl << "Exiting...\n";
    exit(1);
  }
}

double BasisSet::GetNormContractionH ( unsigned int index ) const
{
  if ( index < myNexp )
  {
    return myNContractH.at(index);
  }
  else
  {
    cout << "Invalid index for GetNormContractionH!" << endl << "Exiting...\n";
    exit(1);
  }
}

double BasisSet::GetContractionH ( unsigned int index ) const
{
  if ( index < myNexp )
  {
    return myContractH.at(index);
  }
  else
  {
    cout << "Invalid index for GetContracttionH!" << endl << "Exiting...\n";
    exit(1);
  }
}

double BasisSet::GetExponent(unsigned int index) const
{
  if ( index < myNexp )
  {
    return myExponent.at(index);
  }
  else
  {
    cout << "Invalid index for GetExponent!\nExiting...\n";
    exit(1);
  }
}

unsigned int BasisSet::GetLabel(unsigned int index) const
{
  if ( index < myNBasis )
  {
    return myLabel.at(index);
  }
  else
  {
    cout << "Invalid index for GetLabel!\nExiting...\n";
    exit(1);
  }
}

unsigned int BasisSet::GetPtr(unsigned int index) const
{
  if ( index < myNshell )
  {
    return myPtr.at(index);
  }
  else
  {
    cout << "Invalid index for GetPtr!\nExiting...\n";
    exit(1);
  }
}

void BasisSet::AddSShell ( unsigned int center, const vector<double> & exps,
                           const vector<double> & coeffs )
{
  myNshell++; // Adding a new shell
  myShells.push_back(myNshell-1); // The BF being added is the current shell
  myNcomp.push_back(1); // S shells have one component
  if ( myNshell == 1 ) { myPtr.push_back(0); }
  else { myPtr.push_back(myPtr.at(myNshell-2) + myNprim.at(myNshell-2)); }
  myNprim.push_back(coeffs.size());
  myNexp+=coeffs.size();
  if ( abs( NormContracted(0,coeffs,exps) - 1 ) < 0.001 )
  {
    for ( unsigned int i = 0; i < exps.size(); i++ )
    {
      myExponent.push_back(exps.at(i));
      myNContractS.push_back(coeffs.at(i));
      myNContractP.push_back(0);
      myNContractD.push_back(0);
      myNContractF.push_back(0);
      myNContractG.push_back(0);
      myNContractH.push_back(0);
      myContractS.push_back(coeffs.at(i)/Norm(exps.at(i),1));
      myContractP.push_back(0);
      myContractD.push_back(0);
      myContractF.push_back(0);
      myContractG.push_back(0);
      myContractH.push_back(0);
    }
  }
  else
  {
    for ( unsigned int i = 0; i < exps.size(); i++ )
    {
      myExponent.push_back(exps.at(i));
      myContractS.push_back(coeffs.at(i));
      myContractP.push_back(0);
      myContractD.push_back(0);
      myContractF.push_back(0);
      myContractG.push_back(0);
      myContractH.push_back(0);
      myNContractS.push_back(coeffs.at(i)*Norm(exps.at(i),1));
      myNContractP.push_back(0);
      myNContractD.push_back(0);
      myNContractF.push_back(0);
      myNContractG.push_back(0);
      myNContractH.push_back(0);
    }
  }
}

void BasisSet::AddPShell ( unsigned int center, const vector<double> & exps,
                           const vector<double> & coeffs )
{
  myNshell++; // Adding a new shell
  myShells.push_back(myNshell-1); // The BF being added is the current shell
  myNcomp.push_back(3); // S shells have one component
  if ( myNshell == 1 ) { myPtr.push_back(0); }
  else { myPtr.push_back(myPtr.at(myNshell-2) + myNprim.at(myNshell-2)); }
  myNprim.push_back(coeffs.size());
  myNexp+=coeffs.size();
  if ( abs( NormContracted(1,coeffs,exps) - 1 ) < 0.001 )
  {
    for ( unsigned int i = 0; i < exps.size(); i++ )
    {
      myExponent.push_back(exps.at(i));
      myNContractS.push_back(0);
      myNContractP.push_back(coeffs.at(i));
      myNContractD.push_back(0);
      myNContractF.push_back(0);
      myNContractG.push_back(0);
      myNContractH.push_back(0);
      myContractS.push_back(0);
      myContractP.push_back(coeffs.at(i)/Norm(exps.at(i),103));
      myContractD.push_back(0);
      myContractF.push_back(0);
      myContractG.push_back(0);
      myContractH.push_back(0);
    }
  }
  else
  {
    for ( unsigned int i = 0; i < exps.size(); i++ )
    {
      myExponent.push_back(exps.at(i));
      myNContractS.push_back(0);
      myNContractP.push_back(coeffs.at(i)*Norm(exps.at(i),103));
      myNContractD.push_back(0);
      myNContractF.push_back(0);
      myNContractG.push_back(0);
      myNContractH.push_back(0);
      myContractS.push_back(0);
      myContractP.push_back(coeffs.at(i));
      myContractD.push_back(0);
      myContractF.push_back(0);
      myContractG.push_back(0);
      myContractH.push_back(0);
    }
  }
}

void BasisSet::AddDShell ( unsigned int center, const vector<double> & exps,
                           const vector<double> & coeffs )
{
  myNshell++; // Adding a new shell
  myShells.push_back(myNshell-1); // The BF being added is the current shell
  myNcomp.push_back(5); // S shells have one component
  if ( myNshell == 1 ) { myPtr.push_back(0); }
  else { myPtr.push_back(myPtr.at(myNshell-2) + myNprim.at(myNshell-2)); }
  myNprim.push_back(coeffs.size());
  myNexp+=coeffs.size();
  if ( abs( NormContracted(2,coeffs,exps) - 1 ) < 0.001 )
  {
    for ( unsigned int i = 0; i < exps.size(); i++ )
    {
      myExponent.push_back(exps.at(i));
      myNContractS.push_back(0);
      myNContractP.push_back(0);
      myNContractD.push_back(coeffs.at(i));
      myNContractF.push_back(0);
      myNContractG.push_back(0);
      myNContractH.push_back(0);
      myContractS.push_back(0);
      myContractP.push_back(0);
      myContractD.push_back(coeffs.at(i)/Norm(exps.at(i),255));
      myContractF.push_back(0);
      myContractG.push_back(0);
      myContractH.push_back(0);
    }
  }
  else
  {
    for ( unsigned int i = 0; i < exps.size(); i++ )
    {
      myExponent.push_back(exps.at(i));
      myNContractS.push_back(0);
      myNContractP.push_back(0);
      myNContractD.push_back(coeffs.at(i)*Norm(exps.at(i),255));
      myNContractF.push_back(0);
      myNContractG.push_back(0);
      myNContractH.push_back(0);
      myContractS.push_back(0);
      myContractP.push_back(0);
      myContractD.push_back(coeffs.at(i));
      myContractF.push_back(0);
      myContractG.push_back(0);
      myContractH.push_back(0);
    }
  }
}

void BasisSet::AddFShell ( unsigned int center, const vector<double> & exps,
                           const vector<double> & coeffs )
{
  myNshell++; // Adding a new shell
  myShells.push_back(myNshell-1); // The BF being added is the current shell
  myNcomp.push_back(7); // F shells have 7 components
  if ( myNshell == 1 ) { myPtr.push_back(0); }
  else { myPtr.push_back(myPtr.at(myNshell-2) + myNprim.at(myNshell-2)); }
  myNprim.push_back(coeffs.size());
  myNexp+=coeffs.size();
  if ( abs( NormContracted(3,coeffs,exps) - 1 ) < 0.001 )
  {
    for ( unsigned int i = 0; i < exps.size(); i++ )
    {
      myExponent.push_back(exps.at(i));
      myNContractS.push_back(0);
      myNContractP.push_back(0);
      myNContractD.push_back(0);
      myNContractF.push_back(coeffs.at(i));
      myNContractG.push_back(0);
      myNContractH.push_back(0);
      myContractS.push_back(0);
      myContractP.push_back(0);
      myContractD.push_back(0);
      myContractF.push_back(coeffs.at(i)/Norm(exps.at(i),351));
      myContractG.push_back(0);
      myContractH.push_back(0);
    }
  }
  else
  {
    for ( unsigned int i = 0; i < exps.size(); i++ )
    {
      myExponent.push_back(exps.at(i));
      myNContractS.push_back(0);
      myNContractP.push_back(0);
      myNContractD.push_back(0);
      myNContractF.push_back(coeffs.at(i)*Norm(exps.at(i),351));
      myNContractG.push_back(0);
      myNContractH.push_back(0);
      myContractS.push_back(0);
      myContractP.push_back(0);
      myContractD.push_back(0);
      myContractF.push_back(coeffs.at(i));
      myContractG.push_back(0);
      myContractH.push_back(0);
    }
  }
}

void BasisSet::AddGShell ( unsigned int center, const vector<double> & exps,
                           const vector<double> & coeffs )
{
  myNshell++; // Adding a new shell
  myShells.push_back(myNshell-1); // The BF being added is the current shell
  myNcomp.push_back(9); // G shells have 9 components
  if ( myNshell == 1 ) { myPtr.push_back(0); }
  else { myPtr.push_back(myPtr.at(myNshell-2) + myNprim.at(myNshell-2)); }
  myNprim.push_back(coeffs.size());
  myNexp+=coeffs.size();
  if ( abs( NormContracted(4,coeffs,exps) - 1 ) < 0.001 )
  {
    for ( unsigned int i = 0; i < exps.size(); i++ )
    {
      myExponent.push_back(exps.at(i));
      myNContractS.push_back(0);
      myNContractP.push_back(0);
      myNContractD.push_back(0);
      myNContractF.push_back(0);
      myNContractG.push_back(coeffs.at(i));
      myNContractH.push_back(0);
      myContractS.push_back(0);
      myContractP.push_back(0);
      myContractD.push_back(0);
      myContractF.push_back(0);
      myContractG.push_back(coeffs.at(i)/Norm(exps.at(i),451));
      myContractH.push_back(0);
    }
  }
  else
  {
    for ( unsigned int i = 0; i < exps.size(); i++ )
    {
      myExponent.push_back(exps.at(i));
      myNContractS.push_back(0);
      myNContractP.push_back(0);
      myNContractD.push_back(0);
      myNContractF.push_back(0);
      myNContractG.push_back(coeffs.at(i)*Norm(exps.at(i),451));
      myNContractH.push_back(0);
      myContractS.push_back(0);
      myContractP.push_back(0);
      myContractD.push_back(0);
      myContractF.push_back(0);
      myContractG.push_back(coeffs.at(i));
      myContractH.push_back(0);
    }
  }
}

void BasisSet::AddHShell ( unsigned int center, const vector<double> & exps,
                           const vector<double> & coeffs )
{
  myNshell++; // Adding a new shell
  myShells.push_back(myNshell-1); // The BF being added is the current shell
  myNcomp.push_back(11); // S shells have one component
  if ( myNshell == 1 ) { myPtr.push_back(0); }
  else
  {
    myPtr.push_back(myPtr.at(myNshell-2) + myNprim.at(myNshell-2));
  }
  myNprim.push_back(coeffs.size());
  myNexp+=coeffs.size();
  if ( abs( NormContracted(5,coeffs,exps) - 1 ) < 0.001 )
  {
    for ( unsigned int i = 0; i < exps.size(); i++ )
    {
      myExponent.push_back(exps.at(i));
      myNContractS.push_back(0);
      myNContractP.push_back(0);
      myNContractD.push_back(0);
      myNContractF.push_back(0);
      myNContractG.push_back(0);
      myNContractH.push_back(coeffs.at(i));
      myContractS.push_back(0);
      myContractP.push_back(0);
      myContractD.push_back(0);
      myContractF.push_back(0);
      myContractG.push_back(0);
      myContractH.push_back(coeffs.at(i)/Norm(exps.at(i),551));
    }
  }
  else
  {
    for ( unsigned int i = 0; i < exps.size(); i++ )
    {
      myExponent.push_back(exps.at(i));
      myNContractS.push_back(0);
      myNContractP.push_back(0);
      myNContractD.push_back(0);
      myNContractF.push_back(0);
      myNContractG.push_back(0);
      myNContractH.push_back(coeffs.at(i)*Norm(exps.at(i),551));
      myContractS.push_back(0);
      myContractP.push_back(0);
      myContractD.push_back(0);
      myContractF.push_back(0);
      myContractG.push_back(0);
      myContractH.push_back(coeffs.at(i));
    }
  }
}

void BasisSet::AddBF (unsigned int label, unsigned int center,
                      const vector<double> & exps,
                      const vector<double> & coeffs )
{
  myNBasis++;
  myLabel.push_back(label);
  myCenter.push_back(center);
  // if the label is one of those listed below, we not only add a Basis
  // Function, but we also add a shell
  if ( label == 1 || label == 103 || label ==255 || label == 351 ||
       label == 451 || label == 551 )
  {
    if ( label == 1 ) AddSShell(center,exps,coeffs);
    else if ( label == 103 ) AddPShell(center,exps,coeffs);
    else if ( label == 255 ) AddDShell(center,exps,coeffs);
    else if ( label == 351 ) AddFShell(center,exps,coeffs);
    else if ( label == 451 ) AddGShell(center,exps,coeffs);
    else if ( label == 551 ) AddHShell(center,exps,coeffs);
  }
  else
  {
    // Here we assum that the shells are defined 'together'
    // Seems reasonable that all file formats will obey this
    unsigned int tmp = myShells.at(myShells.size()-1);
    myShells.push_back(tmp);
  }
}

unsigned int BasisSet::GetShell ( unsigned int index ) const
{
  if ( index < myNBasis )
  {
    return myShells.at(index);
  }
  else
  {
    cout << "Invalid index for GetShell!\nExiting...\n";
    exit(1);
  }
}

BasisSet::BasisSet ()
{
  myNshell = 0; myNexp = 0; myNBasis=0;
}

BasisSet::BasisSet ( const BasisSet & obas )
{
  myNBasis = obas.GetNBasis();
  myNshell = obas.GetNShell();
  myNexp = obas.GetNExp();
  myNprim.resize(myNshell);
  myNcomp.resize(myNshell);
  myCenter.resize(myNBasis);
  myPtr.resize(myNshell);
  myLabel.resize(myNBasis);
  myExponent.resize(myNexp);
  myContractS.resize(myNexp);
  myNContractS.resize(myNexp);
  myContractP.resize(myNexp);
  myNContractP.resize(myNexp);
  myContractD.resize(myNexp);
  myNContractD.resize(myNexp);
  myContractF.resize(myNexp);
  myNContractF.resize(myNexp);
  myContractG.resize(myNexp);
  myNContractG.resize(myNexp);
  myContractH.resize(myNexp);
  myNContractH.resize(myNexp);
  for ( unsigned int i = 0; i < myNshell; i++ )
  {
    myNprim.at(i) = obas.GetNumPrim(i);
    myNcomp.at(i) = obas.GetNumComponent(i);
    myPtr.at(i) = obas.GetPtr(i);
  }
  for ( unsigned int i = 0; i < myNBasis; i++ )
  {
    myCenter.at(i) = obas.GetCenter(i);
    myLabel.at(i) = obas.GetLabel(i);
    myShells.at(i) = obas.GetShell(i);
  }
  for ( unsigned int i = 0; i < myNexp; i++ )
  {
    myExponent.at(i) = obas.GetExponent(i);
    myContractS.at(i) = obas.GetContractionS(i);
    myNContractS.at(i) = obas.GetNormContractionS(i);
    myContractP.at(i) = obas.GetContractionP(i);
    myNContractP.at(i) = obas.GetNormContractionP(i);
    myContractD.at(i) = obas.GetContractionD(i);
    myNContractD.at(i) = obas.GetNormContractionD(i);
    myContractF.at(i) = obas.GetContractionF(i);
    myNContractF.at(i) = obas.GetNormContractionF(i);
    myContractG.at(i) = obas.GetContractionG(i);
    myNContractG.at(i) = obas.GetNormContractionG(i);
    myContractH.at(i) = obas.GetContractionH(i);
    myNContractH.at(i) = obas.GetNormContractionH(i);
  }
}

BasisSet & BasisSet::operator=( const BasisSet & obas )
{
  myNBasis = obas.GetNBasis();
  myNshell = obas.GetNShell();
  myNexp = obas.GetNExp();
  myNprim.resize(myNshell);
  myNcomp.resize(myNshell);
  myCenter.resize(myNBasis);
  myPtr.resize(myNshell);
  myLabel.resize(myNBasis);
  myExponent.resize(myNexp);
  myContractS.resize(myNexp);
  myNContractS.resize(myNexp);
  myContractP.resize(myNexp);
  myNContractP.resize(myNexp);
  myContractD.resize(myNexp);
  myNContractD.resize(myNexp);
  myContractF.resize(myNexp);
  myNContractF.resize(myNexp);
  myContractG.resize(myNexp);
  myNContractG.resize(myNexp);
  myContractH.resize(myNexp);
  myNContractH.resize(myNexp);
  for ( unsigned int i = 0; i < myNshell; i++ )
  {
    myNprim.at(i) = obas.GetNumPrim(i);
    myNcomp.at(i) = obas.GetNumComponent(i);
    myPtr.at(i) = obas.GetPtr(i);
  }
  for ( unsigned int i = 0; i < myNBasis; i++ )
  {
    myCenter.at(i) = obas.GetCenter(i);
    myLabel.at(i) = obas.GetLabel(i);
    myShells.at(i) = obas.GetShell(i);
  }
  for ( unsigned int i = 0; i < myNexp; i++ )
  {
    myExponent.at(i) = obas.GetExponent(i);
    myContractS.at(i) = obas.GetContractionS(i);
    myNContractS.at(i) = obas.GetNormContractionS(i);
    myContractP.at(i) = obas.GetContractionP(i);
    myNContractP.at(i) = obas.GetNormContractionP(i);
    myContractD.at(i) = obas.GetContractionD(i);
    myNContractD.at(i) = obas.GetNormContractionD(i);
    myContractF.at(i) = obas.GetContractionF(i);
    myNContractF.at(i) = obas.GetNormContractionF(i);
    myContractG.at(i) = obas.GetContractionG(i);
    myNContractG.at(i) = obas.GetNormContractionG(i);
    myContractH.at(i) = obas.GetContractionH(i);
    myNContractH.at(i) = obas.GetNormContractionH(i);
  }
}
