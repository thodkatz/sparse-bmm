CXX:=g++
CCMPI:=mpic++
CXXFLAGS:=-g -pedantic -Wall

SRC:=src/*.cpp
INC:=inc

$(shell mkdir -p bin)

serial: $(SRC)
	$(CXX) $(CXXFLAGS) -DSERIAL $^ -I$(INC) -o bin/$@

openmp: $(SRC)
	$(CXX) $(CXXFLAGS) -fopenmp -DOPENMP $^ -I$(INC)  -o bin/$@

hybrid: $(SRC)
	$(CCMPI) $(CXXFLAGS) -DHYBRID $^ -I$(INC) -o bin/$@ 

.PHONY: clean

clean:
	rm -f bin/*
