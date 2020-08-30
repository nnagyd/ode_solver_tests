//============================================================================
// Name        : keller_miksis_noroot.cpp
// Author      : Lambert Plavecz
// Version     :
// Copyright   : no
// Description : Parameter study of the Keller-Miksis equation with odeint rkck54
//============================================================================


#include <iostream>
#include <fstream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <string>
#include <chrono>
#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;

const double ro_L = 9.970639504998557e+02;
const double P_inf = 1.0e+5;
const double p_v = 3.166775638952003e+03;
const double sigma = 0.071977583160056;
const double R_E = 10.0/1.0e6;
const double gam = 1.4;
const double c_L = 1.497251785455527e+03;
const double mu_L = 8.902125058209557e-04;
const double theta = 0.0;


int num = 30720;

string file_name = "kellermiksis_cpu_noroot_output.txt";

typedef double value_type;
typedef std::vector< value_type > state_type;

std::vector<double> extrema(2 * num);

class km{
	state_type C;

public:
	km(value_type P_A1, value_type P_A2, value_type omega, value_type omega2){
		C = state_type(13);
		double twr = 2*M_PI/(R_E*omega);
		C[0] = (P_inf - p_v + 2*sigma/R_E)/ro_L*twr*twr;
		C[1] = (1-3*gam)/(ro_L*c_L)*(P_inf - p_v + 2*sigma/R_E)*twr;
		C[2] = (P_inf - p_v)/ro_L * twr*twr;
		C[3] = 2*sigma/(ro_L*R_E) *twr*twr;
		C[4] = 4*mu_L/(ro_L*R_E*R_E) * 2*M_PI/omega;
		C[5] = P_A1/ro_L * twr*twr;
		C[6] = P_A2/ro_L *twr*twr;
		C[7] = R_E * omega*P_A1/(ro_L*c_L) * twr*twr;
		C[8] = R_E * omega*P_A2/(ro_L*c_L) * twr*twr;
		C[9] = R_E*omega/(2*M_PI*c_L);
		C[10] = 3*gam;
		C[11] = omega2/omega;
		C[12] = theta;
	}
	void operator() (const state_type &x, state_type &dxdt, const double t){
		double rx0 = 1.0/x[0];
		double N = (C[0]+C[1]*x[1])*pow(rx0,C[10]) - C[2]*(1.0+C[9]*x[1]) -C[3]*rx0 -C[4]*x[1]*rx0 - (1.0 - C[9]*x[1]/3.0)*1.5*x[1]*x[1]
						-(C[5]*sin(2.0*M_PI*t) + C[6]*sin(2.0*M_PI*C[11]*t + C[12])) * (1.0+C[9]*x[1])
						-x[0]*(C[7]*cos(2.0*M_PI*t) + C[8]*cos(2.0*M_PI*C[11]*t+C[12]));

		double D = x[0] - C[9]*x[0]*x[1] + C[4]*C[9];

		dxdt[0] = x[1];
		dxdt[1] = N/D;
	}
};

class observer
{
    double extr_prev;					 
    state_type x_prev;
    int count = 0;
    int no;

public:
    observer(int row): no(row) { 
    	extrema[2*no] = 0.0; //initial max
		extrema[2*no+1] = 10000.0; //initial min
		extr_prev = 200.0; //some arbitrary large number										  
    }

    void operator()(const state_type &x , double t )
    {
    	if(t==0.0){
    		x_prev = x;
    		return;
    	}
    	if(x[1]*x_prev[1] < 0){//extremum
    		if(count > 2047){ //saving
    			if(x[0] > extrema[2*no]) extrema[2*no] = x[0];
    		   	if(x[0] < extrema[2*no+1]) extrema[2*no+1] = x[0];
				if(count == 2048+63) throw 0; //end integration
			}
    		count++;
			extr_prev = x[0];
    	}
    	x_prev = x;
    }
};
int nums[11] = {256, 768, 1536, 3072, 	3840, 	5120, 	7680, 15360, 30720, 46080, 61440}; // 76800, 92160, 122880, 184320, 307200, 768000, 4147200};
//milliseconds 78032 236444 473084 946772 1184320 1577270 2465706 x 10207031
int main() {
	cout << "Begin" << endl;

	runge_kutta_cash_karp54<state_type> rk_stepper;
	auto stepper = make_controlled(1e-10, 1e-10, rk_stepper) ;
	
	state_type x(2);
	
for(int jj=0; jj < 9;jj++){ //parameter numbers loop
	
	num = nums[jj];
	cout << num << endl;
	extrema = std::vector<double>(2 * num);

	auto t1 = chrono::high_resolution_clock::now();

	double B = 20.0;
	double E = 1000.0;
	double EpB = E/B;
	double invnum = 1.0/(num-1);
	for(int u = 0;u < num;u++){ //f1
		//if(u%100==0) cout << u << endl;
		double f1 = B*pow(EpB, u*invnum);
		x[0] = 1.0;
		x[1] = 0.0;

		km equation( 1.5*1e5, 0.0, 2000.0*M_PI*f1, 0.0);
		observer obs(u);

		try{
			integrate_adaptive(stepper, equation, x, 0.0, 1.0e23, 0.01, obs);
		}catch(...){
			//cout << "Enough" << endl;
		}
	}

	ofstream ofs(file_name);
	if(!ofs.is_open())exit(-1);
	ofs.precision(17);
	//ofs.flags(ios::scientific);
	for(int u = 0;u < num;u++){
		ofs << 1.5 << " " << B*pow(EpB, u*invnum) << " " << 0.0 << " " << 0.0 << " " << theta << " " << R_E
			<< " " << extrema[2*u] << " " << extrema[2*u+1] << "\n";
	}

	ofs.flush();
	ofs.close();

	auto t2 = chrono::high_resolution_clock::now();
	cout << "Done" << endl;
	cout << "Time (ms):" << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << endl;
	
}
	return 0;
}
