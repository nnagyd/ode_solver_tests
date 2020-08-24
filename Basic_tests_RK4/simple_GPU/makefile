CMPL_OPT = -O3 -std=c++11 --ptxas-options=-v --gpu-architecture=sm_35 -lineinfo -maxrregcount=32
SOURCE   = simple_GPU.cu

all: simple_GPU.exe

simple_GPU.exe: $(SOURCE)
	nvcc -o	simple_GPU.exe $(SOURCE) $(CMPL_OPT)

clean:
	rm -f simple_GPU.exe