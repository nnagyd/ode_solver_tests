//============================================================================
// Name        : keller_miksis.cpp
// Author      : Lambert Plavecz
// Version     :
// Copyright   : no
// Description : Parameter study of the Keller-Miksis equation with odeint rkck54
//					includes basic rootfinding algorithm for event handling
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

const double root_tolerance = 1.0e-6;

int num = 3000;

string file_name = "kellermiksis_cpu_output.txt";

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

struct event_exception{
	int event;
	double t_next;
	double t_start;
	state_type x_prev;
	state_type x_next;
};

class observer
{
    double extr_prev;
	double t_prev = 0.0;
	state_type x_prev;
    int count = 0;
    int no;
	bool event = false;

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
		if(event){
			event = false;
			if(count >= 2048){ //saving
    	    	if(x[0] > extrema[2*no]) extrema[2*no] = x[0];
    		   	if(x[0] < extrema[2*no+1]) extrema[2*no+1] = x[0];
    	    	if(count == 2048+63) throw 0; //end integration
    	    }
			count++;
			if(fabs(extr_prev - x[0]) < 1.0e-9 && fabs(x[0]) > 1.0e-9){ //check whether solution converges
				//if true save the current value and end integration
				if(x[0] > extrema[2*no]) extrema[2*no] = x[0]; //max
    			if(x[0] < extrema[2*no+1]) extrema[2*no+1] = x[0]; //min
				throw -1; //end integration
			}
			extr_prev = x[0];  //save extremum for convergence check
			x_prev = x;
			t_prev = t;
			return;
		}
		if(fabs(x[1]) < root_tolerance){
			if(count >= 2048){ //saving
    			if(x[0] > extrema[2*no]) extrema[2*no] = x[0];
    		   	if(x[0] < extrema[2*no+1]) extrema[2*no+1] = x[0];
				if(count == 2048+63) throw 0; //end integration
    		}
    		count++;
    		if(fabs(extr_prev - x[0]) < 1.0e-9 && fabs(x[0]) > 1.0e-9){ //check whether solution converges
    			if(x[0] > extrema[2*no]) extrema[2*no] = x[0]; //max
    			if(x[0] < extrema[2*no+1]) extrema[2*no+1] = x[0]; //min
    	    	throw -1;
    	    }
			extr_prev = x[0];
    	}else if(x[1]*x_prev[1] < 0.0 && fabs(x_prev[1]) > root_tolerance){//extremum
			if(count < 2048){
				count++;
				if(fabs(extr_prev - x[0]) < 1.0e-9 && fabs(x[0]) > 1.0e-9){ //check whether solution converges
					if(x[0] > extrema[2*no]) extrema[2*no] = x[0]; //max
					if(x[0] < extrema[2*no+1]) extrema[2*no+1] = x[0]; //min
					throw -1;
				}
				extr_prev = x[0];
			}else{
				event = true;
				event_exception ie;
				ie.x_prev = x_prev;
				ie.x_next = x;
				ie.t_start = t_prev;
				ie.t_next = t;
				throw ie;
			}
    	}
    	x_prev = x;
		t_prev = t;
    }
};
int nums[11] = {256, 768, 1536, 3072, 	3840, 	5120, 	7680, 15360, 30720, 46080, 61440}; // 76800, 92160, 122880, 184320, 307200, 768000, 4147200};
//milliseconds 59690 179142 358619 717147 896913 1217862 1826911
int main() {
	cout << "Begin" << setprecision(17) << endl;

	typedef runge_kutta_cash_karp54<state_type> stepper_type;
	auto stepper = make_controlled(1.0e-10, 1.0e-10, stepper_type());
	state_type x(2);
	
for(int jj=6; jj < 11;jj++){ //parameter loop
	
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
		
		
		double t_start = 0.0;
		double dt_start = 0.01;
		while(true){
			try{
				integrate_adaptive(stepper, boost::ref(equation), x, t_start, 1.0e23, dt_start, boost::ref(obs));
			}catch(event_exception & ee){
				
				t_start = ee.t_start;
				double t_start_tmp = t_start;
				double t_end = ee.t_next;
				state_type x_start = ee.x_prev;
				state_type x_tmp = x_start;				
				double dt = (t_end-t_start)/2;
				double dt_temp = dt;
				stepper.try_step(boost::ref(equation), x_tmp, t_start_tmp, dt_temp);
				
				while(fabs(x_tmp[1]) > root_tolerance && dt!=0.0){
					if(x_start[1]*x_tmp[1] < 0.0){
						dt = dt/2;
						dt_temp = dt;
						t_start_tmp = t_start;
						x_tmp = x_start;
						stepper.try_step(boost::ref(equation), x_tmp, t_start_tmp, dt_temp);
					}else{
						t_start += dt;
						dt = dt/2;
						dt_temp = dt;
						t_start_tmp = t_start;
						x_start = x_tmp;
						stepper.try_step(boost::ref(equation), x_tmp, t_start_tmp, dt_temp);
					}
				}
				t_start += dt;
				x = x_tmp;
				dt_start = t_end - t_start_tmp;
			}catch(int t){
				//cout << "Enough" << endl;
				break;
			}
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
