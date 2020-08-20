#include <iostream>
#include <ctime>

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

inline void F(double *xp, double x, double p)
{
	*xp = x*x-p;
}

struct vars //4  [double]
{
	double kSum, kAct, x, p;
};

struct integrationConstans //4 [double]
{
	double dt,dtp2,dtp6, c2; //dt = step size, 1/2*dt, 1/6*dt
};

int main()
{
	const int numberOfProblems = 1<<16;

	struct vars v;
	struct integrationConstans consts;

	double x0 = -0.5;
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

	for (int i = 0; i < numberOfProblems; i ++) //parameter sweep loop
	{
		v.p = p_Parameters[i]; //loading parameters from alligned memory
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
			std::cout << "p = " <<v.p[0] << "x =" << v.x << std::endl;
		}

	}//end of parameter sweep

	SimulationEnd = clock() - SimulationStart;
	std::cout << "Elapsed time: " << SimulationEnd/1000 << " ms" << std::endl;
}
