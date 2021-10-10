CXX = g++
CXXFLAGS = -std=c++11 -fopenmp -g -Ofast -funroll-loops # -mavx2 -march=native # edit march depending on system
LDFLAGS = -g -lgomp

BIN = argonmd

CPP_FILES = $(wildcard *.cpp)
OBJ_FILES = $(patsubst %.cpp, %.o, $(CPP_FILES))
.SECONDARY: $(OBJ_FILES)


.PHONY: exe
exe: $(BIN).x

%.x: $(OBJ_FILES)
	$(CXX) $(LDFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<


.PHONY: clean_obj
clean_obj:
	rm -f *.o

.PHONY: clean
clean:
	rm -f *.x *.o
