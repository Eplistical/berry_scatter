CXX = g++ 
MPICXX = mpicxx
OPT = -O2
LIBS += -lboost_program_options -lfftw3 -lopenblas -llapack -lpthread -lgfortran 
SRC = .

all : exact  fssh_conner_mpi fssh_conner_2bath_mpi

exact: $(SRC)/exact.cpp 
	$(CXX) $(OPT) $< -o $@ $(LIBS)

fssh_conner_mpi: $(SRC)/fssh_conner_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

fssh_conner_2bath_mpi: $(SRC)/fssh_conner_2bath_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

check: $(SRC)/check_potential.cpp 
	$(CXX) $(OPT) $< -o $@ $(LIBS)
