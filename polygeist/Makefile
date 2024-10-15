# Usage required CGEIST=cgeist CXX=clang++ CGEIST_LLVM_PATH=Polygeist/llvm-project/install/lib/clang/<version> CGEIST_OPENMP_PATH=Polygeist/llvm-project/openmp/runtime/src MUDOCK=mudock POLY_HEADER=polygeist/include

# Compiler flags
CGEISTFLAGS = -fopenmp --cuda-lower --cpuify=distribute --cuda-gpu-arch=sm_80 -c -std=c++20 --function=* -DEXT_COMPILER -O3 --inner-serialize=1 -resource-dir=$(CGEIST_LLVM_PATH) -DMUDOCK_ENABLE_POLY
CGEISTINCLUDE = -I$(MUDOCK)/include -I${CGEIST_OPENMP_PATH} -I${POLY_HEADER} -DMUDOCK_ENABLE_POLY
CXXFLAGS = -fopenmp -std=c++20 
CXXINCLUDE = -I$(MUDOCK)/include -I${POLY_HEADER}

# Source files for the library
SRC_GPU = src/cuda_batch_sizer.cu \
           src/virtual_screen.cu

SRC_CPU = src/set_device.cpp \
					src/virtual_screen_transpilation.cpp

# Target library name
TARGET = libmudock_poly.a

# Object files
OBJ_GPU = $(SRC_GPU:.cu=.o)
OBJ_CPU = $(SRC_CPU:.cpp=.o)
OBJ = $(OBJ_GPU) $(OBJ_CPU)

# Default target
all: $(TARGET)

# Rule to create the static library
$(TARGET): $(OBJ) src/cuda_object.cu.o src/cuda_object.cpp.o
	llvm-ar rcs $@ $^

# Compile source files into object files
src/cuda_object.cu.o: $(MUDOCK)/src/cuda_implementation/cuda_object.cu
	$(CGEIST) $(CGEISTFLAGS) $(CGEISTINCLUDE) -c $<  -o src/cuda_object.cu.o

# Compile source files into object files
%.o: %.cu
	$(CGEIST) $(CGEISTFLAGS) $(CGEISTINCLUDE) -c $< -o $@

# Compile C++ source files into object files
src/cuda_object.cpp.o: src/cuda_object.cpp
	$(CXX) $(CXXFLAGS) $(CXXINCLUDE) -c $< -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(CXXINCLUDE) -c $< -o $@

# Clean rule
clean:
	rm -f $(OBJ)

clean-all:
	rm -f $(OBJ) $(TARGET)

.PHONY: all clean
