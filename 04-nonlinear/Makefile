CXX=mpic++
CXXFLAGS=-g3 -Wall -O3

HYPRE_SRC_PATH=/home/nico/Desktop/courses/IHPS/mpna/hypre/src

INCLUDE_FLAGS=-I$(HYPRE_SRC_PATH)/hypre/include
LIB_FLAGS=-L$(HYPRE_SRC_PATH)/hypre/lib -lHYPRE

bin=linear_diffusion

compile:
	$(CXX) $(CXXFLAGS) main.cpp $(INCLUDE_FLAGS) $(LIB_FLAGS) -o $(bin)

clean:
	rm $(bin)
