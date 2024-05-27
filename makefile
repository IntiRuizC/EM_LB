CXX=mpic++
CXXFLAGS=-O3
SOURCES = main.cpp Lattice_dynamics.cpp Lattice_Init.cpp vector.cpp
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = Dielectric.x
NP ?= 4

all: $(EXECUTABLE)


$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(OBJECTS) $(LDFLAGS) -o $@


%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@


run: $(EXECUTABLE)
	mpirun -np $(NP) ./$(EXECUTABLE) > pulse.dat


clean:
	rm -f $(OBJECTS) $(EXECUTABLE) *.dat
