## Makefile for INTEGER

ifndef CXX
	CXX := g++
endif
CXXFLAGS += -Wall -O2 -std=c++14
# LDFLAGS += ""

.PHONY: main clean

main: INTEGER

INTEGER: INTEGER.cpp
	$(CXX) $< -o $@ $(CXXFLAGS) $(LDFLAGS)

INTEGER_naive: INTEGER_naive.cpp
	$(CXX) $< -o $@ $(CXXFLAGS) $(LDFLAGS)

INTEGER_map: INTEGER_map.cpp
	$(CXX) $< -o $@ $(CXXFLAGS) $(LDFLAGS)

test: INTEGER
	./$<

clean:
	@rm -f INTEGER INTEGER_naive INTEGER_map
