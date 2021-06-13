CXX = g++
#CXXFLAGS = -O3
CXXFLAGS = -g -O0 -std=c++11 # -Wall -fbounds-check


CPP_FILES = $(wildcard *.cpp)
OBJ_FILES = $(patsubst %.cpp, %.o, $(CPP_FILES))
# .PHONY: variables
# variables:
# 	@echo CPP_FILES: $(CPP_FILES)
# 	@echo OBJ_FILES: $(OBJ_FILES)


.PHONY: all
all: argonmd.x

%.x: $(OBJ_FILES)
	$(CXX) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<


.PHONY: clean
clean:
	rm -f *.x
