CXX = g++ 
MPICXX = mpicxx
OPT = -O3
LIBS += -lboost_program_options -lfftw3 -lopenblas -llapack -lpthread -lgfortran 
SRC = .

all : exact 

exact: $(SRC)/exact.cpp 
	$(CXX) $(OPT) $< -o $@ $(LIBS)
