## Makefile for SinglePoint

CXX := g++
CXXFLAGS += -Wall -O3 -std=c++14
# LDFLAGS  += $(HEJ_LIBS) $(FASTJET_LIBS)

SinglePoint: SinglePoint.cpp
	g++ SinglePoint.cpp $(CXXFLAGS) $(CPPFLAGS) -I. $(LDFLAGS) -o SinglePoint

test: SinglePoint
	./SinglePoint

clean:
	@rm -f SinglePoint
