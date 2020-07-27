#ifndef PERTHREAD_SYSTEMDEFINITION_H
#define PERTHREAD_SYSTEMDEFINITION_H

#define PI 3.14159265358979323846

// SYSTEM
__device__ void PerThread_OdeFunction(int tid, int NT, double* F, double* X, double T, double* cPAR, double* sPAR, int* sPARi, double* ACC, int* ACCi)
{
	double rx1 = 1.0/X[0];
	double p   = pow(rx1, cPAR[8]);
	
	double s1;
	double c1;
	sincospi(2.0*T, &s1, &c1);
	
	double N;
	double D;
	double rD;
	
	N = (cPAR[0]+cPAR[1]*X[1])*p - cPAR[2]*(1.0+cPAR[7]*X[1]) - cPAR[3]*rx1 - cPAR[4]*X[1]*rx1 - 1.5*(1.0-cPAR[7]*X[1]*(1.0/3.0))*X[1]*X[1] - (cPAR[5]*s1) * (1.0+cPAR[7]*X[1]) - X[0]*cPAR[6]*c1;
	D = X[0] - cPAR[7]*X[0]*X[1] + cPAR[4]*cPAR[7];
	rD = 1.0/D;
	
	F[0] = X[1];
	F[1] = N*rD;
}

// EVENTS
__device__ void PerThread_EventFunction(int tid, int NT, double* EF, double* X, double T, double* cPAR, double* sPAR, int* sPARi, double* ACC, int* ACCi)
{
	
}

__device__ void PerThread_ActionAfterEventDetection(int tid, int NT, int IDX, int CNT, double &T, double &dT, double* TD, double* X, double* cPAR, double* sPAR, int* sPARi, double* ACC, int* ACCi)
{
	
}

// ACCESSORIES
__device__ void PerThread_ActionAfterSuccessfulTimeStep(int tid, int NT, double T, double dT, double* TD, double* X, double* cPAR, double* sPAR, int* sPARi, double* ACC, int* ACCi)
{
	if ( X[0]>ACC[0] )
		ACC[0] = X[0];
}

__device__ void PerThread_Initialization(int tid, int NT, double T, double &dT, double* TD, double* X, double* cPAR, double* sPAR, int* sPARi, double* ACC, int* ACCi)
{
	ACC[0] = X[0];
}

__device__ void PerThread_Finalization(int tid, int NT, double T, double dT, double* TD, double* X, double* cPAR, double* sPAR, int* sPARi, double* ACC, int* ACCi)
{
	
}

#endif