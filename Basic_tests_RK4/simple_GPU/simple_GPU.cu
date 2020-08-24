#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <fstream>

using namespace std;

void Linspace(double*, double, double, int);
void Uniform(double*, double, int);

__global__ void RungeKuttaStepOriginal(double*, double*, int);
__device__ void RightHandSide(double&, double, double);

template <class DataType>
DataType* AllocateHostMemory(int);
template <class DataType>
DataType* AllocateDeviceMemory(int);

int main()
{
// INITIAL SETUP ----------------------------------------------------------------------------------
	
	int NumberOfProblems = 61440; // 92160
	int BlockSize        = 128;
	
	cudaSetDevice(1);
	
	double* h_State      = AllocateHostMemory<double>(NumberOfProblems);
	double* h_Parameters = AllocateHostMemory<double>(NumberOfProblems);
	double* d_State      = AllocateDeviceMemory<double>(NumberOfProblems);
	double* d_Parameters = AllocateDeviceMemory<double>(NumberOfProblems);
	
	Linspace(h_Parameters, 0.1, 1.0, NumberOfProblems);
	Uniform(h_State, -0.5, NumberOfProblems);
	
	cudaMemcpy(d_State, h_State, sizeof(double)*NumberOfProblems, cudaMemcpyHostToDevice);
	cudaMemcpy(d_Parameters, h_Parameters, sizeof(double)*NumberOfProblems, cudaMemcpyHostToDevice);
	
	
	int GridSize = NumberOfProblems/BlockSize + (NumberOfProblems % BlockSize == 0 ? 0:1);
	
	clock_t SimulationStart;
	clock_t SimulationEnd;
	
	SimulationStart = clock();
	RungeKuttaStepOriginal<<<GridSize, BlockSize>>> (d_State, d_Parameters, NumberOfProblems);
	SimulationEnd = clock();
	
	cout << "Simulation time: " << 1000.0*(SimulationEnd-SimulationStart) / CLOCKS_PER_SEC << "ms" << endl << endl;
	cout << "Simulation time / 1000 RK4 step: " << 1000.0*(SimulationEnd-SimulationStart) / CLOCKS_PER_SEC << "ms" << endl;
	cout << "Ensemble size:                   " << NumberOfProblems << endl << endl;
		
	cudaMemcpy(h_State, d_State, sizeof(double)*NumberOfProblems, cudaMemcpyDeviceToHost);
	
	//for (int i=0; i<NumberOfProblems; i++)
	//	cout << "P: " << h_Parameters[i] << "   Sates: " << h_State[i] << endl;
}

// AUXILIARY FUNCTION -----------------------------------------------------------------------------

void Linspace(double* x, double B, double E, int N)
{
    double Increment;
	
	x[0]   = B;
	
	if ( N>1 )
	{
		x[N-1] = E;
		Increment = (E-B)/(N-1);
		
		for (int i=1; i<N-1; i++)
		{
			x[i] = B + i*Increment;
		}
	}
}

void Uniform(double* x, double V, int N)
{
	for (int i=0; i<N; i++)
	{
		x[i] = V;
	}
}

__forceinline__ __device__ void RightHandSide(double& F, double X, double P)
{
	F = X*X - P; // 1 FMA
}

__global__ void RungeKuttaStepOriginal(double* d_State, double* d_Parameters, int N)
{
	int tid = threadIdx.x + blockIdx.x*blockDim.x;
	
	if (tid < N)
	{
		double X;
		double P;
		
		double k1;
		double k2;
		double k3;
		double k4;
		double x;
		
		double dT   = 0.01;
		double dTp2 = 0.5*dT;
		double dTp6 = dT * (1.0/6.0);
		
		X = d_State[tid];
		P = d_Parameters[tid];
		
		for (int i=0; i<1000; i++)
		{
			// k1
			RightHandSide(k1, X, P);
			
			x = X + dTp2*k1;
			RightHandSide(k2, x, P);
			
			x = X + dTp2*k2;
			RightHandSide(k3, x, P);
			
			x = X + dT*k3;
			RightHandSide(k4, x, P);
			
			X = X + dTp6*( k1 + 2*k2 + 2*k3 + k4 );
		}
		
		d_State[tid] = X;
	}
}

template <class DataType>
DataType* AllocateHostMemory(int N)
{
    DataType* HostMemory = new (std::nothrow) DataType [N];
    if (HostMemory == NULL)
    {
        std::cerr << "Failed to allocate Memory on the HOST!\n";
        exit(EXIT_FAILURE);
    }
    return HostMemory;
}

template <class DataType>
DataType* AllocateDeviceMemory(int N)
{
    cudaError_t Error = cudaSuccess;
	
	DataType* MemoryAddressInDevice = NULL;
	
	Error = cudaMalloc((void**)&MemoryAddressInDevice, N * sizeof(DataType));
    
	if (Error != cudaSuccess)
    {
        std::cerr << "Failed to allocate Memory on the DEVICE!\n";
        exit(EXIT_FAILURE);
    }
    return MemoryAddressInDevice;
}