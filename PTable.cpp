#include "PTable.h"

using namespace std;

unsigned int PeriodicTable::ElecConfigCount(string el) const
{
  return (SCount(el)+PCount(el)+DCount(el)+FCount(el));
}

unsigned int PeriodicTable::ElecConfigCount(unsigned int atomnum) const
{
  return (SCount(atomnum)+PCount(atomnum)+DCount(atomnum)+FCount(atomnum));
}

unsigned int PeriodicTable::SCount(string el) const
{
  unsigned int anum = GetAtomNumber(el);
  if ( anum < 3 ) return 1;
  else if ( anum < 11 ) return 2;
  else if ( anum < 19 ) return 3;
  else if ( anum < 37 ) return 4;
  else if ( anum < 55 ) return 5;
  else if ( anum < 87 ) return 6;
  else return 7;
}

unsigned int PeriodicTable::FCount(string el) const
{
  unsigned int anum = GetAtomNumber(el);
  if ( anum < 57 ) return 0;
  else if ( anum < 89 ) return 7;
  else return 14;
}

unsigned int PeriodicTable::FCount(unsigned int atomnum) const
{
  if ( atomnum < 57 ) return 0;
  else if ( atomnum < 89 ) return 7;
  else return 14;
}

unsigned int PeriodicTable::DCount(string el) const
{
  unsigned int anum = GetAtomNumber(el);
  if ( anum < 21 ) return 0;
  else if ( anum < 39 ) return 5;
  else if ( anum < 71 ) return 10;
  else if ( anum < 103 ) return 15;
  else return 20;
}

unsigned int PeriodicTable::DCount(unsigned int atomnum) const
{
  if ( atomnum < 21 ) return 0;
  else if ( atomnum < 39 ) return 5;
  else if ( atomnum < 71 ) return 10;
  else if ( atomnum < 103 ) return 15;
  else return 20;
}

unsigned int PeriodicTable::PCount(string el) const
{
  unsigned int anum = GetAtomNumber(el);
  if ( anum < 5 ) return 0;
  else if ( anum < 13 ) return 3;
  else if ( anum < 31 ) return 6;
  else if ( anum < 49 ) return 9;
  else if ( anum < 81 ) return 12;
  else if ( anum < 113 ) return 15;
  else return 18;
}

unsigned int PeriodicTable::PCount(unsigned int atomnum) const
{
  if ( atomnum < 5 ) return 0;
  else if ( atomnum < 13 ) return 3;
  else if ( atomnum < 31 ) return 6;
  else if ( atomnum < 49 ) return 9;
  else if ( atomnum < 81 ) return 12;
  else if ( atomnum < 113 ) return 15;
  else return 18;
}

unsigned int PeriodicTable::SCount(unsigned int atomnum) const
{
  if ( atomnum < 3 ) return 1;
  else if ( atomnum < 11 ) return 2;
  else if ( atomnum < 19 ) return 3;
  else if ( atomnum < 37 ) return 4;
  else if ( atomnum < 55 ) return 5;
  else if ( atomnum < 87 ) return 6;
  else return 7;
}


string PeriodicTable::GetElement(unsigned int atomnum) const
{
  if ( atomnum != 0 )
  {
    return(get<1>(myAtoms.at(atomnum-1)));
  }
  else
  {
    return "X";
  }
}

string PeriodicTable::GetElement(string el) const
{
  if ( el.length() > 2 )
  {
    string upel;
    locale loc;
    stringstream ss;
    for ( string::size_type i = 0; i < el.length(); ++i )
    {
      ss << toupper(el[i],loc);
    }
    ss >> upel;
    for ( unsigned int i = 0; i < myAtoms.size(); i++ )
    {
      if ( get<2>(myAtoms.at(i)) == upel )
      {
        return get<1>(myAtoms.at(i));
      }
    }
  }
  else
  {
    for ( unsigned int i = 0; i < myAtoms.size(); i++ )
    {
      if ( get<1>(myAtoms.at(i)) == el )
      {
        return get<2>(myAtoms.at(i));
      }
    }
  }
}

unsigned int PeriodicTable::GetAtomNumber(string el) const
{
  if ( el.length() > 2 )
  {
    string upel;
    locale loc;
    stringstream ss;
    for ( string::size_type i = 0; i < el.length(); ++i )
    {
      ss << toupper(el[i],loc);
    }
    ss >> upel;
    for ( unsigned int i = 0; i < myAtoms.size(); i++ )
    {
      if ( get<2>(myAtoms.at(i)) == upel )
      {
        return get<0>(myAtoms.at(i));
      }
    }
  }
  else
  {
    for ( unsigned int i = 0; i < myAtoms.size(); i++ )
    {
      if ( get<1>(myAtoms.at(i)) == el )
      {
        return get<0>(myAtoms.at(i));
      }
    }
  }
}

float PeriodicTable::GetWeight(unsigned int atomnum) const
{
  if ( atomnum != 0 )
  {
    return(get<3>(myAtoms.at(atomnum-1)));
  }
  else
  {
    return (get<3>(myAtoms.at(myAtoms.size()-1)));
  }
}

float PeriodicTable::GetWeight(string el) const
{
  if ( el.length() > 2 )
  {
    string upel;
    locale loc;
    stringstream ss;
    for ( string::size_type i = 0; i < el.length(); ++i )
    {
      ss << toupper(el[i],loc);
    }
    ss >> upel;
    for ( unsigned int i = 0; i < myAtoms.size(); i++ )
    {
      if ( get<2>(myAtoms.at(i)) == upel )
      {
        return get<3>(myAtoms.at(i));
      }
    }
  }
  else
  {
    for ( unsigned int i = 0; i < myAtoms.size(); i++ )
    {
      if ( get<1>(myAtoms.at(i)) == el )
      {
        return get<3>(myAtoms.at(i));
      }
    }
  }
}

PeriodicTable::PeriodicTable()
{
  myAtoms.push_back(make_tuple(1,"H","HYDROGEN",1.008));
  myAtoms.push_back(make_tuple(2,"He","HELIUM",4.003));
  myAtoms.push_back(make_tuple(3,"Li","LITHIUM",6.94));
  myAtoms.push_back(make_tuple(4,"Be","BERYLLIUM",9.012));
  myAtoms.push_back(make_tuple(5,"B","BORON",10.81));
  myAtoms.push_back(make_tuple(6,"C","CARBON",12.011));
  myAtoms.push_back(make_tuple(7,"N","NITROGEN",14.007));
  myAtoms.push_back(make_tuple(8,"O","OXYGEN",15.999));
  myAtoms.push_back(make_tuple(9,"F","FLUORINE",18.998));
  myAtoms.push_back(make_tuple(10,"Ne","NEON",20.1798));
  myAtoms.push_back(make_tuple(11,"Na","SODIUM",22.989));
  myAtoms.push_back(make_tuple(12,"Mg","MAGNESIUM",24.305));
  myAtoms.push_back(make_tuple(13,"Al","ALUMINIUM",26.981));
  myAtoms.push_back(make_tuple(14,"Si","SILICON",28.085));
  myAtoms.push_back(make_tuple(15,"P","PHOSPHORUS",30.973));
  myAtoms.push_back(make_tuple(16,"S","SULFUR",32.06));
  myAtoms.push_back(make_tuple(17,"Cl","CHLORINE",35.45));
  myAtoms.push_back(make_tuple(18,"Ar","ARGON",39.948));
  myAtoms.push_back(make_tuple(19,"K","POTASSIUM",39.0983));
  myAtoms.push_back(make_tuple(20,"Ca","CALCIUM",40.078));
  myAtoms.push_back(make_tuple(21,"Sc","SCANDIUM",44.955908));
  myAtoms.push_back(make_tuple(22,"Ti","TITANIUM",47.867));
  myAtoms.push_back(make_tuple(23,"V","VANADIUM",50.9415));
  myAtoms.push_back(make_tuple(24,"Cr","CHROMIUM",51.9961));
  myAtoms.push_back(make_tuple(25,"Mn","MANGANESE",54.938));
  myAtoms.push_back(make_tuple(26,"Fe","IRON",55.845));
  myAtoms.push_back(make_tuple(27,"Co","COBALT",58.933));
  myAtoms.push_back(make_tuple(28,"Ni","NICKEL",58.6934));
  myAtoms.push_back(make_tuple(29,"Cu","COPPER",63.546));
  myAtoms.push_back(make_tuple(30,"Zn","ZINC",65.38));
  myAtoms.push_back(make_tuple(31,"Ga","GALLIUM",69.723));
  myAtoms.push_back(make_tuple(32,"Ge","GERMANIUM",72.63));
  myAtoms.push_back(make_tuple(33,"As","ARSENIC",74.921));
  myAtoms.push_back(make_tuple(34,"Se","SELENIUM",78.971));
  myAtoms.push_back(make_tuple(35,"Br","BROMINE",79.904));
  myAtoms.push_back(make_tuple(36,"Kr","KRYPTON",83.798));
  myAtoms.push_back(make_tuple(37,"Rb","RUBIDIUM",85.4678));
  myAtoms.push_back(make_tuple(38,"Sr","STRONTIUM",87.62));
  myAtoms.push_back(make_tuple(39,"YY","TTRIUM",88.90584));
  myAtoms.push_back(make_tuple(40,"Zr","ZIRCONIUM",91.224));
  myAtoms.push_back(make_tuple(41,"Nb","NIOBIUM",92.906));
  myAtoms.push_back(make_tuple(42,"Mo","MOLYBDENUM",95.95));
  myAtoms.push_back(make_tuple(43,"Tc","TECHNETIUM",97));
  myAtoms.push_back(make_tuple(44,"Ru","RUTHENIUM",101.07));
  myAtoms.push_back(make_tuple(45,"Rh","RHODIUM",102.90549));
  myAtoms.push_back(make_tuple(46,"Pd","PALLADIUM",106.42));
  myAtoms.push_back(make_tuple(47,"Ag","SILVER",107.8682));
  myAtoms.push_back(make_tuple(48,"Cd","CADMIUM",112.414));
  myAtoms.push_back(make_tuple(49,"In","INDIUM",114.818));
  myAtoms.push_back(make_tuple(50,"Sn","TIN",118.711));
  myAtoms.push_back(make_tuple(51,"Sb","ANTIMONY",121.76));
  myAtoms.push_back(make_tuple(52,"Te","TELLURIUM",127.6));
  myAtoms.push_back(make_tuple(53,"I","IODINE",126.904));
  myAtoms.push_back(make_tuple(54,"Xe","XENON",131.293));
  myAtoms.push_back(make_tuple(55,"Cs","CAESIUM",132.905451));
  myAtoms.push_back(make_tuple(56,"Ba","BARIUM",137.327));
  myAtoms.push_back(make_tuple(57,"La","LANTHANUM",138.90547));
  myAtoms.push_back(make_tuple(58,"Ce","CERIUM",140.116));
  myAtoms.push_back(make_tuple(59,"Pr","PRASEODYMIUM",140.90766));
  myAtoms.push_back(make_tuple(60,"Nd","NEODYMIUM",144.242));
  myAtoms.push_back(make_tuple(61,"Pm","PROMETHIUM",145));
  myAtoms.push_back(make_tuple(62,"Sm","SAMARIUM",150.36));
  myAtoms.push_back(make_tuple(63,"Eu","EUROPIUM",151.964));
  myAtoms.push_back(make_tuple(64,"Gd","GADOLINIUM",157.25));
  myAtoms.push_back(make_tuple(65,"Tb","TERBIUM",158.925355));
  myAtoms.push_back(make_tuple(66,"Dy","DYSPROSIUM",162.5));
  myAtoms.push_back(make_tuple(67,"Ho","HOLMIUM",164.93));
  myAtoms.push_back(make_tuple(68,"Er","ERBIUM",167.259));
  myAtoms.push_back(make_tuple(69,"Tm","THULIUM",168.934219));
  myAtoms.push_back(make_tuple(70,"Yb","YTTERBIUM",173.045));
  myAtoms.push_back(make_tuple(71,"Lu","LUTETIUM",174.9668));
  myAtoms.push_back(make_tuple(72,"Hf","HAFNIUM",178.49));
  myAtoms.push_back(make_tuple(73,"Ta","TANTALUM",180.947));
  myAtoms.push_back(make_tuple(74,"W","TUNGSTEN",183.84));
  myAtoms.push_back(make_tuple(75,"Re","RHENIUM",186.207));
  myAtoms.push_back(make_tuple(76,"Os","OSMIUM",190.23));
  myAtoms.push_back(make_tuple(77,"Ir","IRIDIUM",192.217));
  myAtoms.push_back(make_tuple(78,"Pt","PLATINUM",195.085));
  myAtoms.push_back(make_tuple(79,"Au","GOLD",196.966));
  myAtoms.push_back(make_tuple(80,"Hg","MERCURY",200.592));
  myAtoms.push_back(make_tuple(81,"Tl","THALLIUM",204.38));
  myAtoms.push_back(make_tuple(82,"Pb","LEAD",207.2));
  myAtoms.push_back(make_tuple(83,"Bi","BISMUTH",208.98));
  myAtoms.push_back(make_tuple(84,"Po","POLONIUM",209));
  myAtoms.push_back(make_tuple(85,"At","ASTATINE",210));
  myAtoms.push_back(make_tuple(86,"Rn","RADON",222));
  myAtoms.push_back(make_tuple(87,"Fr","FRANCIUM",223));
  myAtoms.push_back(make_tuple(88,"Ra","RADIUM",226));
  myAtoms.push_back(make_tuple(89,"Ac","ACTINIUM",227));
  myAtoms.push_back(make_tuple(90,"Th","THORIUM",232.0377));
  myAtoms.push_back(make_tuple(91,"Pa","PROTACTINIUM",231.035));
  myAtoms.push_back(make_tuple(92,"U","URANIUM",238.028));
  myAtoms.push_back(make_tuple(93,"Np","NEPTUNIUM",237));
  myAtoms.push_back(make_tuple(94,"Pu","PLUTONIUM",244));
  myAtoms.push_back(make_tuple(95,"Am","AMERICIUM",243));
  myAtoms.push_back(make_tuple(96,"Cm","CURIUM",247));
  myAtoms.push_back(make_tuple(97,"Bk","BERKELIUM",247));
  myAtoms.push_back(make_tuple(98,"Cf","CALIFORNIUM",251));
  myAtoms.push_back(make_tuple(99,"Es","EINSTEINIUM",252));
  myAtoms.push_back(make_tuple(100,"Fm","FERMIUM",257));
  myAtoms.push_back(make_tuple(101,"Md","MENDELEVIUM",258));
  myAtoms.push_back(make_tuple(102,"No","NOBELIUM",259));
  myAtoms.push_back(make_tuple(103,"Lr","LAWRENCIUM",262));
  myAtoms.push_back(make_tuple(0,"X","DUMMY",-1));
}
