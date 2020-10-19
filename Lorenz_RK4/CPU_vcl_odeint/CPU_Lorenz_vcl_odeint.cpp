#include <array>
#include <chrono>
#include <boost/numeric/odeint.hpp>
#include "vectorclass.h"

using namespace std;
using namespace boost::numeric::odeint;

const double num = 40320;
const double inv_num = 1/num;
const int step_number = 1000;
const double dt = 0.01;
const double t_max = dt*step_number;
const int unroll = 16;
const int vectorSize = 4;
const int arraySize = 3*unroll;

typedef std::array<Vec4d, arraySize> state_type;

const double sigma = 10;
const double beta = 8.0/3.0;

Vec4d k_global[unroll];


template<int unroll>
void lorenz(const state_type &x, state_type &dxdt, double t)
{
	for(int i = 0, j = 0; i < unroll; i++, j+=3)
	{
		dxdt[j] = sigma*(x[j+1] - x[j+0]);
		dxdt[j+1] = x[j]*(k_global[i] - x[j+2]) - x[j+1];
		dxdt[j+2] = x[j]*x[j+1] - beta*x[j+2];
	}
}

int main() {
	runge_kutta4_classic<state_type> stepper;
	state_type x;

	auto t1 = chrono::high_resolution_clock::now();

	for(int u = 0;u < num;u+=vectorSize*unroll)
	{
		for(int ii=0;ii < arraySize;ii++)
		{
			x[ii] = 10;
		}

		for(int ii = 0,j = 00; ii < vectorSize*unroll; ii+=vectorSize,j++)
		{
			k_global[j] = Vec4d((u+ii) *   21.0 * inv_num,
								(u+ii+1) * 21.0 * inv_num,
								(u+ii+2) * 21.0 * inv_num,
								(u+ii+3) * 21.0 * inv_num);
		}

		integrate_const(stepper, lorenz<unroll>, x, 0.0, t_max, dt);
	}

	auto t2 = chrono::high_resolution_clock::now();
	cout << "Done" << endl;
	cout << "Time (ms):" << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << endl;
	return 0;
}
