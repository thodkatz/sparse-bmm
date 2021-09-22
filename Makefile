CXX:=g++
CCMPI:=mpic++
CXXFLAGS:=-O3 -pedantic -Wall

SRC:=src/*.cpp

$(shell mkdir -p bin)

serial: $(SRC)
	$(CXX) $(CXXFLAGS) -DSERIAL $^ -o bin/$@

mpi: $(SRC)
	$(MPICC) $(CXXFLAGS) -DMPI $^ -o bin/$@

hybrid: $(SRC)
	$(MPICC) $(CXXFLAGS) -fopenmp -DHYBRID $^ -o bin/$@ 

.PHONY: clean

clean:
	rm -f bin/*