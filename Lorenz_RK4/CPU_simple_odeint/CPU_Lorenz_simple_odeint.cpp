//============================================================================
// Name        : lorenz_RK4.cpp
// Author      : Lambert Plavecz
// Version     :
// Copyright   : no
// Description : lorenz with odeint RK4
//============================================================================

//rollout 8 is fastest
#include <array>
#include <chrono>
#include <boost/numeric/odeint.hpp>
#include <boost/convert.hpp>

using namespace std;
using namespace boost::numeric::odeint;

int num = 32720;
double inv_num = 1.0/(num-1);
const int step_number = 1000;
const double dt = 0.01;
const double t_max = dt*step_number;
const int rollOut = 8;
const int arraySize = 3 * rollOut;

typedef std::array< double, arraySize> state_type;

const double sigma = 10;
const double beta = 8.0/3.0;

double k_global = 0.0;
double k_global2 = 0.0;
double k_global3 = 0.0;
double k_global4 = 0.0;
double k_global5 = 0.0;
double k_global6 = 0.0;
double k_global7 = 0.0;
double k_global8 = 0.0;
void lorenz(const state_type &x, state_type &dxdt, double t){
	dxdt[0] = sigma*(x[1] - x[0]);
	dxdt[1] = x[0]*(k_global - x[2]) - x[1];
	dxdt[2] = x[0]*x[1] - beta*x[2];
	dxdt[3] = sigma*(x[1+3] - x[0+3]);
	dxdt[4] = x[0+3]*(k_global2 - x[2+3]) - x[1+3];
	dxdt[5] = x[0+3]*x[1+3] - beta*x[2+3];
	dxdt[6] = sigma*(x[1+6] - x[0+6]);
	dxdt[7] = x[0+6]*(k_global2 - x[2+6]) - x[1+6];
	dxdt[8] = x[0+6]*x[1+6] - beta*x[2+6];
	dxdt[9] = sigma*(x[1+9] - x[0+9]);
	dxdt[10] = x[0+9]*(k_global2 - x[2+9]) - x[1+9];
	dxdt[11] = x[0+9]*x[1+9] - beta*x[2+9];
	dxdt[12] = sigma*(x[1+12] - x[0+12]);
	dxdt[13] = x[0+12]*(k_global2 - x[2+12]) - x[1+12];
	dxdt[14] = x[0+12]*x[1+12] - beta*x[2+12];
	dxdt[15] = sigma*(x[1+15] - x[0+15]);
	dxdt[16] = x[0+15]*(k_global2 - x[2+15]) - x[1+15];
	dxdt[17] = x[0+15]*x[1+15] - beta*x[2+15];
	dxdt[18] = sigma*(x[1+18] - x[0+18]);
	dxdt[19] = x[0+18]*(k_global2 - x[2+18]) - x[1+18];
	dxdt[20] = x[0+18]*x[1+18] - beta*x[2+18];
	dxdt[21] = sigma*(x[1+21] - x[0+21]);
	dxdt[22] = x[0+21]*(k_global2 - x[2+21]) - x[1+21];
	dxdt[23] = x[0+21]*x[1+21] - beta*x[2+21];
}

int nums[13] = {256, 512,  1024,  2048,  4096, 7680, 15360, 46080, 92160, 184320, 307200, 768000, 4147200};
int main() {
	cout << "Lorenz RK4 started" << endl;
	runge_kutta4<state_type, double, state_type, double, array_algebra> stepper;
	state_type x;
	
	/*for(int j = 0;j < 13;j++){
		num = nums[j];
		inv_num = 1.0/(num-1);*/
		cout << num << endl;
		auto t1 = chrono::high_resolution_clock::now();
	
		for(int u = 0;u < num;u+=rollOut){
			for(int ii=0;ii < arraySize;ii++) x[ii] = 10;

			k_global = u * 21 * inv_num;
			k_global2 = (u+1) * 21 * inv_num;
			k_global3 = (u+2) * 21 * inv_num;
			k_global4 = (u+3) * 21 * inv_num;
			k_global5 = (u+4) * 21 * inv_num;
			k_global6 = (u+5) * 21 * inv_num;
			k_global7 = (u+6) * 21 * inv_num;
			k_global8 = (u+7) * 21 * inv_num;

			integrate_const(stepper, lorenz, x, 0.0, t_max, dt);

			//cout << "x0: " << x[0] << ", x1: " << x[1] << ", x2: " << x[2] << endl;
		}

		auto t2 = chrono::high_resolution_clock::now();
		cout << "Done" << endl;
		cout << "Time (ms):" << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << endl;
	//}
	
	return 0;
}
