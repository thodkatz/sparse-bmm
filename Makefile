CXX:=g++
CCMPI:=mpic++
CXXFLAGS:=-g -pedantic -Wall

SRC:=src/*.cpp
INC:=inc

$(shell mkdir -p bin)

serial: $(SRC)
	$(CXX) $(CXXFLAGS) -DSERIAL $^ -I$(INC) -o bin/$@

mpi: $(SRC)
	$(MPICC) $(CXXFLAGS) -DMPI $^ -o bin/$@

hybrid: $(SRC)
	$(MPICC) $(CXXFLAGS) -fopenmp -DHYBRID $^ -o bin/$@ 

.PHONY: clean

clean:
	rm -f bin/*