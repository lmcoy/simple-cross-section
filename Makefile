CXX=g++
CXXFLAGS=-std=c++11 -O2 -Wall `lhapdf-config --cppflags`
LDFLAGS=`lhapdf-config --ldflags`

OBJS=vegas.o fourmomentum.o rnd.o parameters_sm.o uux_mupmum.o libaloha.a phasespace.o

.PHONY: all clean

all: xsec ${OBJS}

libaloha.a: 
	cd aloha; make; cp libaloha.a ..

vegas.o: vegas_mpi.c vegas.h
	mpicc -c $< -o $@

fourmomentum.o: fourmomentum.cpp fourmomentum.h
	${CXX} ${CXXFLAGS} -c -o $@ $<

rnd.o: rnd.cpp rnd.h
	${CXX} ${CXXFLAGS} -c -o $@ $<

parameters_sm.o: parameters_sm.cpp parameters_sm.h
	${CXX} ${CXXFLAGS} -c -o $@ $<

uux_mupmum.o: uux_mupmum.cpp uux_mupmum.h parameters_sm.o
	${CXX} ${CXXFLAGS} -c -o $@ $<

phasespace.o: phasespace.cpp phasespace.h
	${CXX} ${CXXFLAGS} -c -o $@ $<

xsec: main.cpp ${OBJS}
	mpicxx ${LDFLAGS} ${CXXFLAGS} -o $@ $^

clean:
	rm -rf ${OBJS} xsec