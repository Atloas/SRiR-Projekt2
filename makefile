NN=2

all: nbody run clean

nbody: nbody.cpp
	UPCXX_GASNET_CONDUIT=udp upcxx -O nbody.cpp -o nbody

run: nbody
	upcxx-run -n $(NN) $$(upcxx-nodes nodes) nbody

.PHONY: all run clean setup

clean:
	rm -f nbody resultdata.txt