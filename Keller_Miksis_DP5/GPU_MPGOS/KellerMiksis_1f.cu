#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>

#include "KellerMiksis_1f.cuh"
#include "SingleSystem_PerThread.cuh"

#define PI 3.14159265358979323846

using namespace std;

// Physical control parameters
const int NumberOfFrequency = 128;
	
// Solver Configuration
#define SOLVER RKCK45 // RK4, RKCK45
const int NT   = NumberOfFrequency; // NumberOfThreads
const int SD   = 2;     // SystemDimension
const int NCP  = 9;    // NumberOfControlParameters
const int NSP  = 0;     // NumberOfSharedParameters
const int NISP = 0;     // NumberOfIntegerSharedParameters
const int NE   = 0;     // NumberOfEvents
const int NA   = 1;     // NumberOfAccessories
const int NIA  = 0;     // NumberOfIntegerAccessories
const int NDO  = 0;     // NumberOfPointsOfDenseOutput

void Linspace(vector<double>&, double, double, int);
void Logspace(vector<double>&, double, double, int);
void FillSolverObject(ProblemSolver<NT,SD,NCP,NSP,NISP,NE,NA,NIA,NDO,SOLVER,double>&, const vector<double>&, int);

int main()
{
	int BlockSize = 64;
	
	vector<double> Frequency(NT,0);
	Logspace(Frequency, 20.0, 1000.0, NT);
	
	// Setup CUDA a device
	ListCUDADevices();
	
	int MajorRevision  = 3;
	int MinorRevision  = 5;
	int SelectedDevice = SelectDeviceByClosestRevision(MajorRevision, MinorRevision);
	
	PrintPropertiesOfSpecificDevice(SelectedDevice);
	
	// Setup Solver
	ProblemSolver<NT,SD,NCP,NSP,NISP,NE,NA,NIA,NDO,SOLVER,double> ScanKellerMiksis(SelectedDevice);
	
	ScanKellerMiksis.SolverOption(ThreadsPerBlock, BlockSize);
	ScanKellerMiksis.SolverOption(RelativeTolerance, 0, 1e-10);
	ScanKellerMiksis.SolverOption(RelativeTolerance, 1, 1e-10);
	ScanKellerMiksis.SolverOption(AbsoluteTolerance, 0, 1e-10);
	ScanKellerMiksis.SolverOption(AbsoluteTolerance, 1, 1e-10);
	
// SIMULATIONS ------------------------------------------------------------------------------------
	
	vector<double> GlobalMaxima(NT,0);
	FillSolverObject(ScanKellerMiksis, Frequency, NT);
	
	clock_t SimulationStart = clock();
	ScanKellerMiksis.SynchroniseFromHostToDevice(All);
	for (int i=0; i<1024; i++)
	{
		ScanKellerMiksis.Solve();
		ScanKellerMiksis.InsertSynchronisationPoint();
		ScanKellerMiksis.SynchroniseSolver();
	}
	ScanKellerMiksis.SynchroniseFromDeviceToHost(All);
	ScanKellerMiksis.InsertSynchronisationPoint();
	ScanKellerMiksis.SynchroniseSolver();
	clock_t TransientSimulationEnd = clock();
		cout << "Transient simulation time: " << (double)(TransientSimulationEnd-SimulationStart) / CLOCKS_PER_SEC << "s" << endl << endl;
	
	for (int i=0; i<64; i++)
	{
		ScanKellerMiksis.Solve();
		ScanKellerMiksis.SynchroniseFromDeviceToHost(Accessories);
		ScanKellerMiksis.InsertSynchronisationPoint();
		ScanKellerMiksis.SynchroniseSolver();
		
		for (int tid=0; tid<NT; tid++)
			GlobalMaxima[tid] = fmax( ScanKellerMiksis.GetHost(tid, Accessories, 0), GlobalMaxima[tid] );
	}
	
	// Save collected data to file
	ofstream DataFile;
	DataFile.open ( "KellerMiksis_1f.txt" );
	int Width = 18;
	DataFile.precision(10);
	DataFile.flags(ios::scientific);
	
	for (int tid=0; tid<NT; tid++)
	{
		DataFile.width(8); DataFile << tid << ',';
		DataFile.width(Width); DataFile << Frequency[tid] << ',';
		DataFile.width(Width); DataFile << ScanKellerMiksis.GetHost(tid, ActualState, 0) << ',';
		DataFile.width(Width); DataFile << ScanKellerMiksis.GetHost(tid, ActualState, 1) << ',';
		DataFile.width(Width); DataFile << GlobalMaxima[tid];
		DataFile << '\n';
	}
	
	DataFile.close();
	clock_t TotalSimulationEnd = clock();
		cout << "Total simulation time: " << (double)(TotalSimulationEnd-SimulationStart) / CLOCKS_PER_SEC << "s" << endl << endl;
	
	//ScanKellerMiksis.Print(DenseOutput, 25);
}

// ------------------------------------------------------------------------------------------------

void Linspace(vector<double>& x, double B, double E, int N)
{
    double Increment;
	
	x[0] = B;
	
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

void Logspace(vector<double>& x, double B, double E, int N)
{
    x[0] = B; 
	
	if ( N>1 )
	{
		x[N-1] = E;
		double ExpB = log10(B);
		double ExpE = log10(E);
		double ExpIncr = (ExpE-ExpB)/(N-1);
		for (int i=1; i<N-1; i++)
		{
			x[i] = pow(10,ExpB + i*ExpIncr);
		}
	}
}

// ------------------------------------------------------------------------------------------------

void FillSolverObject(ProblemSolver<NT,SD,NCP,NSP,NISP,NE,NA,NIA,NDO,SOLVER,double>& Solver, const vector<double>& Values, int NumberOfThreads)
{
	// Declaration of physical control parameters
	double P2; // frequency          [kHz]
	
	// Declaration of constant parameters
	double P1 =  1.5; // pressure amplitude   [bar]
	double P6 = 10.0; // equilibrium radius   [mum]
	double P7 =  1.0; // ambient pressure     [bar]
	double P9 =  1.4; // polytrophic exponent [-]
	
	// Material properties
	double Pv  = 3.166775638952003e+03;
    double Rho = 9.970639504998557e+02;
    double ST  = 0.071977583160056;
    double Vis = 8.902125058209557e-04;
    double CL  = 1.497251785455527e+03;
	
	// Auxiliary variables
	double Pinf;
	double PA1;
	double RE;
	double f1;
	
	for (int i=0; i<NumberOfThreads; i++)
	{	
		// Update physical parameters
		P2 = Values[i]; // frequency [kHz]
		
		Solver.SetHost(i, TimeDomain, 0, 0);
		Solver.SetHost(i, TimeDomain, 1, 1);
		
		// Initial conditions are the equilibrium condition y1=1; y2=0;
		Solver.SetHost(i, ActualState, 0, 1.0);
		Solver.SetHost(i, ActualState, 1, 0.0);
		
		// Scaling of physical parameters to SI
		Pinf = P7 * 1e5;
		PA1  = P1 * 1e5;
		RE   = P6 / 1e6;
		
		// Scale to angular frequency
		f1   = 2.0*PI*(P2*1000);
					
		// System coefficients and other, auxiliary parameters
		Solver.SetHost(i, ControlParameters, 0, (2.0*ST/RE + Pinf - Pv) * pow(2.0*PI/RE/f1, 2.0) / Rho );
		Solver.SetHost(i, ControlParameters, 1, (1.0-3.0*P9) * (2*ST/RE + Pinf - Pv) * (2.0*PI/RE/f1) / CL/Rho );
		Solver.SetHost(i, ControlParameters, 2, (Pinf - Pv) * pow(2.0*PI/RE/f1, 2.0) / Rho );
		Solver.SetHost(i, ControlParameters, 3, (2.0*ST/RE/Rho) * pow(2.0*PI/RE/f1, 2.0) );
		Solver.SetHost(i, ControlParameters, 4, (4.0*Vis/Rho/pow(RE,2.0)) * (2.0*PI/f1) );
		Solver.SetHost(i, ControlParameters, 5, PA1 * pow(2.0*PI/RE/f1, 2.0) / Rho );
		Solver.SetHost(i, ControlParameters, 6, (RE*f1*PA1/Rho/CL) * pow(2.0*PI/RE/f1, 2.0) );
		Solver.SetHost(i, ControlParameters, 7, RE*f1/(2.0*PI)/CL );
		Solver.SetHost(i, ControlParameters, 8, 3.0*P9 );
	}
}