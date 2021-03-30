CXX = g++
CXXFLAGS = -O3

.PHONY: all
all: argonmd.x

%.x: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

.PHONY: clean
clean:
	rm -f *.x
	