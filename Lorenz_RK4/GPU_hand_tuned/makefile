INCL_DIR = -I$(HOME)/ode_solver_tests/Lorenz_RK4/GPU_hand_tuned/SourceCodes
CMPL_OPT = -O3 -std=c++11 --ptxas-options=-v --gpu-architecture=sm_35 -lineinfo -maxrregcount=64
SOURCE   = Lorenz_HandTuned.cu

all: Lorenz_HandTuned.exe

Lorenz_HandTuned.exe: $(SOURCE)
	nvcc -o	Lorenz_HandTuned.exe $(SOURCE) $(INCL_DIR) $(CMPL_OPT)

clean:
	rm -f Lorenz_HandTuned.exe