#include "mautil.h"

unsigned int DFac (unsigned int arg)
{
  unsigned int total = 1;
  if ( (arg % 2) == 0 )
  {
    for ( unsigned int i = 1; i <= arg/2; i++ ) total*=(2*i);
  }
  else
  {
    for ( unsigned int i = 1; i<= (arg+1)/2; i++ ) total*=(2*i-1);
  }
  return total;
}

void GoTo ( ifstream & inFile, string searchterm )
{
  string line;
  getline(inFile,line);
  while ( line.find(searchterm) == string::npos )
  {
    getline(inFile,line);
  }
}

unsigned int ShellToLabel ( string shell )
{
  if ( shell == "s" || shell == "S" ) return 1;
  if ( shell == "p" || shell == "P" ) return 103;
  if ( shell == "d" || shell == "D" ) return 255;
  if ( shell == "f" || shell == "F" ) return 351;
  if ( shell == "g" || shell == "G" ) return 451;
  if ( shell == "h" || shell == "H" ) return 551; 
}
