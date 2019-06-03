## Makefile for INTEGER

ifndef CXX
	CXX := g++
endif
ifndef PREFIX
	PREFIX := /usr/bin
endif
ifndef FASTJET_CXX
	FASTJET_CXX := $(shell fastjet-config --cxxflags)
endif
ifndef FASTJET_LIBS
	FASTJET_LIBS := $(shell fastjet-config --libs)
endif

CXXFLAGS += -Wall -O2 -std=c++14 $(FASTJET_CXX)
LDFLAGS += $(FASTJET_LIBS)

OTHER_PROGRAMS = INTEGER_naive INTEGER_map
PROGRAMS = INTEGER $(OTHER_PROGRAMS)

.PHONY: main clean all

main: INTEGER

all: $(PROGRAMS)

INTEGER: %: %.cpp %.hh Cuts.hh
	$(CXX) $< -o $@ $(CXXFLAGS) $(LDFLAGS)

$(OTHER_PROGRAMS): %: %.cpp
	$(CXX) $< -o $@ $(CXXFLAGS) $(LDFLAGS)

test: INTEGER
	./$<

install: INTEGER
	cp $< $(PREFIX)

clean:
	@rm -f $(PROGRAMS)
