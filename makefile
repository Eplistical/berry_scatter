CXX = g++ 
MPICXX = mpicxx
OPT = -O2
LIBS += -lboost_program_options -lfftw3 -lopenblas -llapack -lpthread -lgfortran 
SRC = .

all : exact  fssh

exact: $(SRC)/exact.cpp 
	$(CXX) $(OPT) $< -o $@ $(LIBS)

fssh: $(SRC)/fssh_conner_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

check: $(SRC)/check_potential.cpp 
	$(CXX) $(OPT) $< -o $@ $(LIBS)
