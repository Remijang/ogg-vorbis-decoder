NVFLAGS := -std=c++11 -O3 -use_fast_math -ftz=true
TARGET := main

.PHONY: all
all: $(TARGET)

$(TARGET): mdct.o main.cu
	nvcc $(NVFLAGS) -o main main.cu mdct.o

mdct.o: mdct.cu
	nvcc $(NVFLAGS) -c mdct.cu

clean:
	rm -rf main *.o
