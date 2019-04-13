## Makefile for SinglePoint

ifndef CXX
	CXX := g++
endif
CXXFLAGS += -Wall -O2 -std=c++14
# LDFLAGS += ""

SinglePoint: SinglePoint.cpp
	$(CXX) $< -o $@ $(CXXFLAGS) -I.  $(LDFLAGS)

test: SinglePoint
	./$<

clean:
	@rm -f SinglePoint
