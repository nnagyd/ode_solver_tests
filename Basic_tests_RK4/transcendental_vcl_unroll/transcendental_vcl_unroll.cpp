#include <iostream>
#include <ctime>
#define MAX_VECTOR_SIZE 256
#include "vectorclass.h"
#include "vectormath_trig.h"

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
	*xp = p*sin(x);
}

template<int unroll>
struct vars //16 * unroll [double]
{
	Vec4d kSum[unroll], kAct[unroll], x[unroll], p[unroll];
};

struct integrationConstans //16 [double]
{
	Vec4d dt,dtp2,dtp6,c2; //dt = step size, 1/2*dt, 1/6*dt
};

int main()
{
	const int numberOfProblems = 1<<17;
	const int unroll = 2; //amount of manual unrolling of outer for loop
	const int outerForLoopStep = 4 * unroll;

	struct vars<unroll> v;
	struct integrationConstans consts;

	Vec4d x0 = 0.1;
	double * p_Parameters = linspace(0.1, 1.0, numberOfProblems);
	const int numberOfSteps = 1000;
	const double dt = 0.01;
	consts.dtp6 = dt / 6.;
	consts.dtp2 = dt / 2.;
	consts.dt = dt;

	clock_t SimulationStart;
	clock_t SimulationEnd;
	SimulationStart = clock();

	for (int i = 0; i < numberOfProblems; i += outerForLoopStep) //parameter sweep loop
	{
		for (int j = 0, offset = i; j < unroll; j++, offset += 4) //initial setup before integration
		{
			v.p[j].load(p_Parameters + offset); //loading parameters from alligned memory
			v.x[j] = x0;
		}

		for (int j = 0; j < numberOfSteps; j++) //integration loop
		{
			for (int k = 0; k < unroll; k++) //k1
			{
				F(&(v.kAct[k]), v.x[k], v.p[k]); //k1 = F(x)
				v.kSum[k] = v.kAct[k];
				v.kAct[k] = v.x[k] + consts.dtp2*v.kAct[k];
			}

			for (int k = 0; k < unroll; k++) //k2
			{
				F(&(v.kAct[k]), v.kAct[k], v.p[k]); //k2 = F(x+dt*k1/2)
				v.kSum[k] += consts.c2*v.kAct[k];
				v.kAct[k] = v.x[k] + consts.dtp2*v.kAct[k];
			}

			for (int k = 0; k < unroll; k++) //k3
			{
				F(&(v.kAct[k]), v.kAct[k], v.p[k]); //k3 = F(x+dt*k2/2)
				v.kSum[k] += consts.c2*v.kAct[k];
				v.kAct[k] = v.x[k] + consts.dt*v.kAct[k];
			}

			for (int k = 0; k < unroll; k++) //k4
			{
				F(&(v.kAct[k]), v.kAct[k], v.p[k]); //k4 = F(x+dt*k3)
				v.kSum[k] += v.kAct[k];
			}

			for (int k = 0; k < unroll; k++) //step
			{
				v.x[k] += v.kSum[k]*consts.dtp6; //x = x + xp*dt
			}

		}//end of integration loop

		if (i % 5000 == 1) //check end values
		{
			for (int k = 0; k < unroll; k++)
			{
				std::cout << "p = " << ((double*)(&v.p[k]))[0] << "x =" << ((double*)(&v.x[k]))[0] << std::endl;
			}
		}

	}//end of parameter sweep

	SimulationEnd = clock() - SimulationStart;
	std::cout << "Elapsed ido: " << SimulationEnd/1000 << " ms" << std::endl;
}
