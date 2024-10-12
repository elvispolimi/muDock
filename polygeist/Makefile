CGEIST=/home/gianmarco.accordi/muDock/cgeist
CGEIST_LLVM_PATH=/home/gianmarco.accordi/Polygeist/llvm-project/install/18
CGEIST_OPENMP_PATH=/home/gianmarco.accordi/Polygeist/llvm-project/openmp/runtime/src
MUDOCK_HEADER=$(shell pwd)/../mudock/include

# Compiler flags
CXXFLAGS =-fopenmp --cuda-lower --cpuify=distribute --cuda-gpu-arch=sm_80 -c -std=c++20 --function=* -DEXT_COMPILER -O3 --inner-serialize=1 -resource-dir=$(CGEIST_LLVM_PATH)
INCLUDE=-I$(MUDOCK_HEADER) -I${CGEIST_OPENMP_PATH}

# Source files for the library
SRC_GPU = src/cuda_implementation/cuda_batch_sizer.cu \
           src/cuda_implementation/cuda_object.cu \
           src/cuda_implementation/virtual_screen.cu

SRC_CPU = src/cuda_implementation/set_device.cpp

# Target library name
TARGET = libmudock_poly.a

# Object files
OBJ_GPU = $(SRC_GPU:.cu=.o)
OBJ_CPU = $(SRC_CPU:.cpp=.o)
OBJ = $(OBJ_GPU) $(OBJ_CPU)

# Default target
all: $(TARGET)

# Rule to create the static library
$(TARGET): $(OBJ)
	llvm-ar rcs $@ $^

# Compile source files into object files
%.o: %.cu
	$(CGEIST) $(CXXFLAGS) $(INCLUDE) -c $< -o $@

# Compile C++ source files into object files
%.o: %.cpp
	$(CXX) -fopenmp -I$(MUDOCK_HEADER) -c $< -o $@

# Clean rule
clean:
	rm -f $(OBJ)

.PHONY: all clean
