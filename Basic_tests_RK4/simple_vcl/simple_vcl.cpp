#include <iostream>
#include <ctime>
#define MAX_VECTOR_SIZE 256
#include "vectorclass.h"

double * linspace(double a, double b, int numberOfInts)
{
	double * list = (double*)aligned_alloc(64, numberOfInts * sizeof(double));
	double d = (b - a) / (double(numberOfInts) - 1);
	list[0] = a;
	for (int i = 1; i < numberOfInts; i++)
	{
		list[i] = list[i - 1] + d;
	}
	return list;
}

inline void F(Vec4d *xp, Vec4d x, Vec4d p)
{
	*xp = x*x-p;
}

struct vars //16  [double]
{
	Vec4d kSum, kAct, x, p;
};

struct integrationConstans //16 [double]
{
	Vec4d dt,dtp2,dtp6, c2; //dt = step size, 1/2*dt, 1/6*dt
};

int main()
{
	const int numberOfProblems = 1<<16;
	const int outerForLoopStep = 4;

	struct vars v;
	struct integrationConstans consts;

	Vec4d x0 = -0.5;
	double * p_Parameters = linspace(0.1, 1.0, numberOfProblems);
	const int numberOfSteps = 1000;
	const double dt = 0.01;
	consts.dtp6 = dt / 6.;
	consts.dtp2 = dt / 2.;
	consts.dt = dt;
	consts.c2 = 2.0;

	clock_t SimulationStart;
	clock_t SimulationEnd;
	SimulationStart = clock();

	for (int i = 0; i < numberOfProblems; i += outerForLoopStep) //parameter sweep loop
	{
		v.p.load_a(p_Parameters + i); //loading parameters from alligned memory
		v.x = x0;

		for (int j = 0; j < numberOfSteps; j++) //integration loop
		{
			F(&(v.kAct), v.x, v.p); //k1 = F(x)
			v.kSum = v.kAct;
			v.kAct = v.x + consts.dtp2*v.kAct;

			F(&(v.kAct), v.kAct, v.p); //k2 = F(x+dt*k1/2)
			v.kSum += consts.c2*v.kAct;
			v.kAct = v.x + consts.dtp2*v.kAct;

			F(&(v.kAct), v.kAct, v.p); //k3 = F(x+dt*k2/2)
			v.kSum += consts.c2*v.kAct;
			v.kAct = v.x + consts.dt*v.kAct;

			F(&(v.kAct), v.kAct, v.p); //k4 = F(x+dt*k3)
			v.kSum += v.kAct;

			v.x += v.kSum*consts.dtp6; //x = x + xp*dt
		}//end of integration loop

		if (i % 5000 == 1) //check end values
		{
			std::cout << "p = " << ((double*)(&v.p))[0] << "x =" << ((double*)(&v.x))[0] << std::endl;
		}

	}//end of parameter sweep

	SimulationEnd = clock() - SimulationStart;
	std::cout << "Elapsed time: " << SimulationEnd/1000 << " ms" << std::endl;
}
