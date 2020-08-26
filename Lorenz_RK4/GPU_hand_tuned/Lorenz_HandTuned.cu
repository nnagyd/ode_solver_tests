#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <fstream>

using namespace std;

void Linspace(double*, double, double, int);
void Uniform(double*, double, int);

__global__ void RungeKuttaStepOriginal(double* __restrict__, const double* __restrict__, int);
__global__ void RungeKuttaStepRegisterFriendly(double* __restrict__, const double* __restrict__, int);
__device__ void Lorenz(double* __restrict__, const double* __restrict__, double);

template <class DataType>
DataType* AllocateHostMemory(int);
template <class DataType>
DataType* AllocateDeviceMemory(int);

int main()
{
// INITIAL SETUP ----------------------------------------------------------------------------------
	
	int NumberOfProblems = 768000;
	int NumberOfThreads  = NumberOfProblems;
	int BlockSize        = 64;
	
	cudaSetDevice(1);
	
	double* h_State      = AllocateHostMemory<double>( 3*NumberOfProblems );
	double* h_Parameters = AllocateHostMemory<double>(   NumberOfProblems );
	double* d_State      = AllocateDeviceMemory<double>( 3*NumberOfProblems );
	double* d_Parameters = AllocateDeviceMemory<double>(   NumberOfProblems );
	
	Linspace(h_Parameters, 0.0, 21.0, NumberOfProblems);
	Uniform(h_State, 10.0, NumberOfProblems);
	Uniform(&h_State[   NumberOfProblems ], 10.0, NumberOfProblems);
	Uniform(&h_State[ 2*NumberOfProblems ], 10.0, NumberOfProblems);
	
	
	cudaMemcpy(d_State, h_State, 3*sizeof(double)*NumberOfProblems, cudaMemcpyHostToDevice);
	cudaMemcpy(d_Parameters, h_Parameters, sizeof(double)*NumberOfProblems, cudaMemcpyHostToDevice);
	
	int GridSize = NumberOfThreads/BlockSize + (NumberOfThreads % BlockSize == 0 ? 0:1);
	
	clock_t SimulationStart;
	clock_t SimulationEnd;
	
	SimulationStart = clock();
	RungeKuttaStepRegisterFriendly<<<GridSize, BlockSize>>> (d_State, d_Parameters, NumberOfProblems);
	cudaDeviceSynchronize();
	SimulationEnd = clock();
	
	cout << "Simulation time: " << 1000.0*(SimulationEnd-SimulationStart) / CLOCKS_PER_SEC << "ms" << endl << endl;
	cout << "Simulation time / 1000 RK4 step: " << 1000.0*(SimulationEnd-SimulationStart) / CLOCKS_PER_SEC << "ms" << endl;
	cout << "Ensemble size:                   " << NumberOfProblems << endl << endl;
		
	cudaMemcpyAsync(h_State, d_State, 3*sizeof(double)*NumberOfProblems, cudaMemcpyDeviceToHost);
	
	
	//for (int i=0; i<NumberOfProblems; i++)
	//	cout << "P: " << h_Parameters[i] << "   Sates: " << h_State[i] << ", " << h_State[i+NumberOfProblems] << ", " << h_State[i+2*NumberOfProblems] << endl;
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

__forceinline__ __device__ void Lorenz(double* __restrict__ F, const double* __restrict__ X, double P)
{
	// How 5 FMA and 1 ADD/MUL is possible
	F[0] = 10.0*(X[1] - X[0]); // 2 FP inst: 1 FMA, 1 ADD
	F[1] = P*X[0] - X[1] - X[0]*X[2]; // 2 FP inst: 2 FMA
	F[2] = X[0]*X[1] - 2.666 * X[2]; // 2 FP inst: 1 MUL, 1 FMA
}

__global__ void RungeKuttaStepRegisterFriendly(double* __restrict__ d_State, const double* __restrict__ d_Parameters, int N)
{
	int tid = threadIdx.x + blockIdx.x*blockDim.x;
	
	if (tid < N)
	{
		double X[3];
		double P;
		
		double k1[3];
		double ks[3];
		double x[3];
		
		double T    = 0.0;
		double dT   = 1e-3;
		double dTp2 = 0.5*dT;
		double dTp6 = dT * (1.0/6.0);
		//double t;
		
		X[0] = d_State[tid];
		X[1] = d_State[tid + N];
		X[2] = d_State[tid + 2*N];
		
		P = d_Parameters[tid];
		
		for (int i=0; i<1000; i++)
		{
			// k1
			Lorenz(k1, X, P);
			
			// k2
			//t = T + 0.5*dT;
			
			#pragma unroll 3
			for (int j=0; j<3; j++)
			{
				x[j]  = X[j] + dTp2*k1[j];
				ks[j] = k1[j];
			}
			
			Lorenz(k1, x, P);
			
			// k3
			//t = T + 0.5*dT;
			
			#pragma unroll 3
			for (int j=0; j<3; j++)
			{
				x[j]  = X[j] + dTp2*k1[j];
				ks[j] = ks[j]+2*k1[j];
			}
			
			Lorenz(k1, x, P);
			
			// k4
			//t = T + dT;
			
			#pragma unroll 3
			for (int j=0; j<3; j++)
			{
				x[j] = X[j] + dT*k1[j];
				ks[j] = ks[j]+2*k1[j];
			}
			
			Lorenz(k1, x, P);
			
			
			// Update state
			#pragma unroll 3
			for (int j=0; j<3; j++)
				X[j] = X[j] + dTp6*( ks[j] + k1[j] );
			
			T += dT;
		}
		
		d_State[tid] = X[0];
		d_State[tid + N] = X[1];
		d_State[tid + 2*N] = X[2];
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