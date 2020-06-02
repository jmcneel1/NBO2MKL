TARGET = nbo2mkl

SRC = nbo2mkl.cpp BasisSet.cpp Atom.cpp mautil.cpp PTable.cpp
OBJ = $(SRC:%.cpp=%.o)

INCLUDE = -I

LIB =

CXX = g++
CXXFLAGS = -Wall -g --std=c++11

CALL = $(CXX) $(CXXFLAGS)

all: $(TARGET)

$(TARGET): $(OBJ) $(TARGET).cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(TARGET).cpp $(OBJ) $(LIB) -o $(TARGET)

clean:
	rm -f $(OBJ) $(TARGET).o $(TARGET)

depend:
	makedepend -- $(CXXFLAGS) -- $(INCLUDE) -- $(SRC) $(TARGET).cpp

# ESSENTIAL LINE
nbo2mkl.o: nbo2mkl.cpp BasisSet.h BasisSet.o Atom.h Atom.o mautil.h mautil.o
	$(CALL) -c nbo2mkl.cpp
BasisSet.o: BasisSet.cpp mautil.o mautil.h
	$(CALL) -c BasisSet.cpp
Atom.o: Atom.cpp PTable.o PTable.h
	$(CALL) -c Atom.cpp
PTable.o: PTable.cpp PTable.h
	$(CALL) -c PTable.cpp
mautil.o: mautil.cpp mautil.h
	$(CALL) -c mautil.cpp
