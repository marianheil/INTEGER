## Makefile for INTEGER

ifndef CXX
	CXX := g++
endif
CXXFLAGS += -Wall -O2 -std=c++14
# LDFLAGS += ""

PROGRAMS = INTEGER INTEGER_naive INTEGER_map

.PHONY: main clean all

main: INTEGER

all: $(PROGRAMS)

$(PROGRAMS): %: %.cpp
	$(CXX) $< -o $@ $(CXXFLAGS) $(LDFLAGS)

test: INTEGER
	./$<

clean:
	@rm -f $(PROGRAMS)
