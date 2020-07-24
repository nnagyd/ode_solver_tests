#include <iostream>
#include <ctime>
#define MAX_VECTOR_SIZE 256
#include "vectorclass.h"

double * linspace(double a, double b, int numberOfInts)
{
	double * list = (double*)aligned_alloc(64, numberOfInts * sizeof(double));
	double d = (b - a) / double(numberOfInts - 1);
	list[0] = a;
	for (int i = 1; i < numberOfInts; i++)
	{
		list[i] = list[i - 1] + d;
	}
	return list;
}

template<int unroll>
inline void F(Vec4d *xp, Vec4d *x, Vec4d *p)
{
	for (int i = 0, j = 0; j < unroll; i += 3,j++)
	{
		xp[i] = 10.*(x[i+1] - x[i]);
		xp[i+1] = p[j] * x[i] - x[i+1] - x[i] * x[i+2];
		xp[i+2] = x[i] * x[i+1] - (8./3.) * x[i+2];
	}
}

template<int totalNumberOfVariables,int numberOfEquations>
struct RK4dynamicVariables //4*(4*totalNumberOfVars + numberOfEquations) [double]
{
	Vec4d kAct[totalNumberOfVariables];
	Vec4d kSum[totalNumberOfVariables];
	Vec4d x[totalNumberOfVariables];
	Vec4d xTmp[totalNumberOfVariables];
	Vec4d p[numberOfEquations];
};


struct RK4constants //4*4 [double]
{
	Vec4d c2 = Vec4d(2.);
	Vec4d dt, dtp2, dtp6;
};

int main()
{
	const int numberOfProblems = 46080;
	const int unroll = 2; //amount of manual unrolling of outer for loop
	const int numberOfVariablesPerEquation = 3;
	const int numberOfVariablesTotal = numberOfVariablesPerEquation * unroll;
	const int numberOfSteps = 1000;
	const int outerForLoopStep = 4 * unroll;
	const double dt = 1e-2;
	RK4dynamicVariables<numberOfVariablesTotal, unroll> v;
	RK4constants c;
	c.dt = dt;
	c.dtp2 = dt / 2.;
	c.dtp6 = dt / 6.;

	Vec4d x0[numberOfVariablesTotal];
	for (int i = 0; i < numberOfVariablesTotal; i+=3)
	{
		x0[i] = 10.;
		x0[i+1] = 10.;
		x0[i+2] = 10.;
	}

	double * p_Parameters = linspace(0, 21.0, numberOfProblems);

	clock_t SimulationStart;
	clock_t SimulationEnd;
	SimulationStart = clock();

	for (int i = 0; i < numberOfProblems; i+=outerForLoopStep) //parameter sweep loop
	{
		for (int l = 0,offset = i; l < unroll; l++,offset += 4)
		{
			v.p[l].load_a(p_Parameters + offset); //loading parameters from alligned memory
		}

		for (int l = 0; l < numberOfVariablesTotal; l++)
		{
			v.x[l] = x0[l]; //initial condition
		}

		for (int j = 0; j < numberOfSteps; j++) //integration loop
		{
			F<unroll>(v.kAct, v.x, v.p); //k1
			for (int l = 0; l < numberOfVariablesTotal; l++)
			{
				v.kSum[l] = v.kAct[l];
				v.xTmp[l] = v.x[l] + c.dtp2 * v.kAct[l];
			}

			F<unroll>(v.kAct, v.xTmp, v.p); //k2
			for (int l = 0; l < numberOfVariablesTotal; l++)
			{
				v.kSum[l] += c.c2*v.kAct[l];
				v.xTmp[l] = v.x[l] + c.dtp2 * v.kAct[l];
			}

			F<unroll>(v.kAct, v.xTmp, v.p); //k3
			for (int l = 0; l < numberOfVariablesTotal; l++)
			{
				v.kSum[l] += c.c2*v.kAct[l];
				v.xTmp[l] = v.x[l] + c.dt * v.kAct[l];
			}

			F<unroll>(v.kAct, v.xTmp, v.p); //k4
			for (int l = 0; l < numberOfVariablesTotal; l++)
			{
				v.kSum[l] += v.kAct[l];
				v.x[l] += c.dtp6*v.kSum[l]; //step
			}
		} //end of integration

		if (i % 1000 == 1)
		{
			std::cout << "x1 = " << (v.x[0])[0] << "\tx2 = " << (v.x[1])[0] << "\tx3 = " << (v.x[2])[0] << "\tp= " << (v.p[0])[0] << std::endl;
		}

	}//end of parameter sweep

	SimulationEnd = clock() - SimulationStart;
	std::cout << "Elapsed time: " << SimulationEnd/1000 << " ms" << std::endl;
}
