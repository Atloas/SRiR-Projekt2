NN=2

nbody: nbody.cpp
	UPCXX_GASNET_CONDUIT=udp upcxx -O nbody.cpp -o nbody

all: run

run: nbody
	upcxx-run -n $(NN) $$(upcxx-nodes nodes) nbody

.PHONY: all run
.PHONY: clean

clean:
	rm -f nbody resultdata.txt
