CXX = g++
CXXFLAGS = -std=c++11 -g -O3 # -Wall -fbounds-check

BIN = argonmd

CPP_FILES = $(wildcard *.cpp)
OBJ_FILES = $(patsubst %.cpp, %.o, $(CPP_FILES))


.PHONY: exe
exe: $(BIN).x

%.x: $(OBJ_FILES)
	$(CXX) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<


.PHONY: clean
clean:
	rm -f *.x
