#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <fstream>

#include "SingleSystem_PerThread_IndexingMacroEnabled.cuh"
#include "Lorenz_SystemDefinition.cuh"
#include "SingleSystem_PerThread_IndexingMacroDisabled.cuh"
#include "SingleSystem_PerThread.cuh"

#define PI 3.14159265358979323846

#define SOLVER RK4
#define EVNT   EVNT0
#define DOUT   DOUT0

using namespace std;

void Linspace(vector<double>&, double, double, int);
void FillSolverObject(ProblemSolver<SOLVER,EVNT,DOUT>&, const vector<double>&, int);

int main()
{
// INITIAL SETUP ----------------------------------------------------------------------------------
	
	int NumberOfProblems = 92160;
	int NumberOfThreads  = NumberOfProblems;
	int BlockSize        = 64;
	
	ListCUDADevices();
	
	int MajorRevision  = 3;
	int MinorRevision  = 5;
	int SelectedDevice = SelectDeviceByClosestRevision(MajorRevision, MinorRevision);
	
	PrintPropertiesOfSpecificDevice(SelectedDevice);
	
	int NumberOfParameters_R = NumberOfProblems;
	double R_RangeLower = 0.0;
    double R_RangeUpper = 56.0;
		vector<double> Parameters_R_Values(NumberOfParameters_R,0);
		Linspace(Parameters_R_Values, R_RangeLower, R_RangeUpper, NumberOfParameters_R);
	
	
	ConstructorConfiguration ConfigurationDuffing;
	
	ConfigurationDuffing.NumberOfThreads = NumberOfThreads;
	ConfigurationDuffing.SystemDimension = 3;
	ConfigurationDuffing.NumberOfControlParameters = 1;
	
	ProblemSolver<SOLVER,EVNT,DOUT> ScanDuffing(ConfigurationDuffing, SelectedDevice);
	
	ScanDuffing.SolverOption(ThreadsPerBlock, BlockSize);
	ScanDuffing.SolverOption(InitialTimeStep, 0.001);
	ScanDuffing.SolverOption(MaximumNumberOfTimeSteps, 10000);
	
	
// SIMULATIONS ------------------------------------------------------------------------------------
	
	clock_t SimulationStart;
	clock_t SimulationEnd;
	
	FillSolverObject(ScanDuffing, Parameters_R_Values, NumberOfThreads);
	
	ScanDuffing.SynchroniseFromHostToDevice(All);
	ScanDuffing.InsertSynchronisationPoint();
	ScanDuffing.SynchroniseSolver();
		
	SimulationStart = clock();
		ScanDuffing.Solve();
		ScanDuffing.InsertSynchronisationPoint();
		ScanDuffing.SynchroniseSolver();
	SimulationEnd = clock();
		cout << "Simulation time: " << 1000.0*(SimulationEnd-SimulationStart) / CLOCKS_PER_SEC << "ms" << endl << endl;
	
	ScanDuffing.SynchroniseFromDeviceToHost(All);
	ScanDuffing.InsertSynchronisationPoint();
	ScanDuffing.SynchroniseSolver();
	
	//for (int i=0; i<NumberOfProblems; i++)
	//	cout << ScanDuffing.GetHost(i, ActualState, 0) << endl;
	
	
	cout << "Test finished!" << endl;
}

// AUXILIARY FUNCTION -----------------------------------------------------------------------------

void Linspace(vector<double>& x, double B, double E, int N)
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

void FillSolverObject(ProblemSolver<SOLVER,EVNT,DOUT>& Solver, const vector<double>& R_Values, int NumberOfThreads)
{
	int ProblemNumber = 0;
	for (int k=0; k<NumberOfThreads; k++)
	{
		Solver.SetHost(ProblemNumber, TimeDomain,  0, 0 );
		Solver.SetHost(ProblemNumber, TimeDomain,  1, 0.001*10000 );
		
		Solver.SetHost(ProblemNumber, ActualState, 0, 10.0 );
		Solver.SetHost(ProblemNumber, ActualState, 1, 10.0 );
		Solver.SetHost(ProblemNumber, ActualState, 2, 10.0 );
		
		Solver.SetHost(ProblemNumber, ControlParameters, 0, R_Values[k] );
		
		ProblemNumber++;
	}
}