CC=mpic++
R=mpiexec
NN=4

nbody: nbody.cpp
	$(CC) -o nbody nbody.cpp

all: run

run: nbody
	$(R) -f nodes -n $(NN) ./nbody 

.PHONY: all run
.PHONY: clean

clean:
	rm -f nbody resultdata.txt
