//============================================================================
// Name        : valve.cpp
// Author      : Lambert Plavecz
// Version     : 2
// Copyright   : no
// Description : Parameter study of a pressure relief valve equation with odeint and rkck54 
// 				 includes basic root finding algorithm
//============================================================================


#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <chrono>
#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;

const double kappa = 1.25;
const double beta = 20.0;
const double delta = 10.0;
const double r = 0.8;

const double root_tolerance = 1.0e-6;

const int num = 15360; //1024 - 56s, 30720 - 1201s
const double invnum = 1.0 / (num -1);

string file_name = "impact_dyn_output_root.txt";

double matrix[num][64];

typedef std::vector< double > state_type;

double q = 0.0;

class impact_dyn {
public:
    impact_dyn(){}
    void operator()(const state_type &x, state_type &dxdt, const double t){
    	dxdt[0] = x[1]; //x[0] displacement, x[1] velocity, x[2] pressure in chamber
    	dxdt[1] = -kappa*x[1] - (x[0] + delta) + x[2]; //
    	dxdt[2] = beta * (q - x[0]* sqrt(x[2]));
    }
};

struct event_exception{
	int event;
	double t_next;
	double t_start;
	state_type x_prev;
	state_type x_next;
	
};

class impact_observer
{
    double extr_prev;
	double t_prev = 0.0;
	state_type x_prev;
	int row_number;
    int count = 0;
    int extremum_count = 0;
	bool event = false;
	int event_number = 0;

public:

    impact_observer(int u): row_number(u){
		extr_prev = 200.0; //some arbitrary large number
	}

    void operator()(state_type &x , double t )
    {
    	if(t==0.0){
    		x_prev = x;
    		return;
    	}
		
		if(event){
			
			if(event_number == 0) x[1] = -r * x[1];
			event = false;
			if(count >= 2048){ //saving
    	    	matrix[row_number][extremum_count++] = x[0];
    	    	if( extremum_count == 64) throw 0; //end integration
    	    }
			count++;
			if(fabs(extr_prev - x[0]) < 1.0e-9 && fabs(x[0]) > 1.0e-9){ //check whether solution converges
				//if true save the current value and end integration
				while(extremum_count < 64){ //write out the convergent solution into the array
				    matrix[row_number][extremum_count++] = x[0];
				}
				throw -1; //end integration
			}
			extr_prev = x[0];  //save extremum for convergence check
			x_prev = x;
			t_prev = t;
			return;
		}
		if(fabs(x[0]) < root_tolerance){
			
			x[1] = -r * x[1];
			if(count >= 2048){ //saving
    	    	matrix[row_number][extremum_count++] = x[0];
    	    	if( extremum_count == 64) throw 0; //end integration
    	    }
			count++;
			/*if(fabs(extr_prev - x[0]) < 1.0e-9 && fabs(x[0]) > 1.0e-9){ //check whether solution converges
				//if true save the current value and end integration
				while(extremum_count < 64){ //write out the convergent solution into the array
				    matrix[row_number][extremum_count++] = x[0];
				}
				throw -1; //end integration
			}*/
			extr_prev = x[0];  //save extremum for convergence check
		}else if(fabs(x[1]) < root_tolerance){
			
			if(count >= 2048){ //saving
    	    	matrix[row_number][extremum_count++] = x[0];
    	    	if( extremum_count == 64) throw 0; //end integration
    	    }
			count++;
			if(fabs(extr_prev - x[0]) < 1.0e-9 && fabs(x[0]) > 1.0e-9){ //check whether solution converges
				//if true save the current value and end integration
				while(extremum_count < 64){ //write out the convergent solution into the array
				    matrix[row_number][extremum_count++] = x[0];
				}
				throw -1; //end integration
			}
			extr_prev = x[0];  //save extremum for convergence check
		}else if(x[0]*x_prev[0] < 0.0 && x[1] < 0.0){ //impact
			event = true;
			event_number = 0;
			event_exception ie;
			ie.x_prev = x_prev;
			ie.x_next = x;
			ie.t_start = t_prev;
			ie.t_next = t;
			ie.event = 0;
			throw ie;
    	}
		else if(x[1]*x_prev[1] < 0.0 && fabs(x_prev[1]) > root_tolerance){//extremum
    		event = true;
			event_number = 1;
			event_exception ie;
			ie.x_prev = x_prev;
			ie.x_next = x;
			ie.t_start = t_prev;
			ie.t_next = t;
			ie.event = 1;
			throw ie;
    	}
    	x_prev = x;
		t_prev = t;
    }
};


//param number 256, 768, 1536, 3072, 3840, 5120, 7680, 15360, 30720, 46080, 61440}; // 76800, 92160, 122880, 184320, 307200, 768000, 4147200};
//milliseconds 9618 28925 57665 114692 143240 191204 286491 574039 1201000
int main() {
	double tolerance = 10e-10;
	cout << "Impact dynamics started\n" << setprecision(17) << num << endl;

	typedef runge_kutta_cash_karp54< state_type , double , state_type , double > stepper_type;
	
	state_type x(3);

	auto t1 = chrono::high_resolution_clock::now();

	impact_dyn impact;
	auto stepper = make_controlled( tolerance , tolerance, stepper_type() );
	
	for(int u = 0;u < num;u++){
		//if(u%1000==0)cout << u << endl;
		x[0] = 0.2;
		x[1] = 0.0;
		x[2] = 0.0;

		q = 0.2 + u * 9.8 * invnum;
		impact_observer obs(u);
		double t_start = 0.0;
		while(true){
			try{
				integrate_adaptive(stepper, impact, x, t_start, 1.0e23, 0.01, boost::ref(obs));
			}catch(event_exception & ee){
				
				t_start = ee.t_start;
				double t_start_tmp = t_start;
				double t_end = ee.t_next;
				state_type x_start = ee.x_prev;
				state_type x_tmp = x_start;				
				double dt = (t_end-t_start)/2;
				double dt_temp = dt;
				stepper.try_step(impact, x_tmp, t_start_tmp, dt_temp);
				
				if(ee.event == 0) //impact
				while(fabs(x_tmp[0]) > root_tolerance && dt!=0.0){
					if(x_start[0]*x_tmp[0] < 0.0){
						dt = dt/2;
						dt_temp = dt;
						t_start_tmp = t_start;
						x_tmp = x_start;
						stepper.try_step(impact, x_tmp, t_start_tmp, dt_temp);
					}else{
						t_start += dt;
						dt = dt/2;
						dt_temp = dt;
						t_start_tmp = t_start;
						x_start = x_tmp;
						stepper.try_step(impact, x_tmp, t_start_tmp, dt_temp);
					}
					
				}
				else if(ee.event == 1) //extremum
				while(fabs(x_tmp[1]) > root_tolerance && dt!=0.0){
					if(x_start[1]*x_tmp[1] < 0.0){
						dt = dt/2;
						dt_temp = dt;
						t_start_tmp = t_start;
						x_tmp = x_start;
						stepper.try_step(impact, x_tmp, t_start_tmp, dt_temp);
						//integrate_adaptive(stepper, impact, x_tmp, t_start, t_start+dt, dt/2);
					}else{
						t_start += dt;
						dt = dt/2;
						dt_temp = dt;
						t_start_tmp = t_start;
						x_start = x_tmp;
						stepper.try_step(impact, x_tmp, t_start_tmp, dt_temp);
						//integrate_adaptive(stepper, impact, x_tmp, t_start, t_start+dt, dt/2);
					}
					
				}
				else cout << "Bad event!" << endl;
				t_start += dt;
				x = x_tmp;
			}catch(int t){
				//cout << "Enough" << endl;
				break;
			}
		}
	}

	ofstream ofs(file_name);
	if(!ofs.is_open())exit(-1);
	ofs.precision(17);
	ofs.flags(ios::scientific);
	for(int u = 0;u < num;u++){ //f1
		ofs << 0.2 + u * 9.8 * invnum;
		for(int i=0; i < 64;i++){
			ofs << " " << matrix[u][i];
		}
		ofs << "\n";
	}

	ofs.flush();
	ofs.close();
	
	auto t2 = chrono::high_resolution_clock::now();
	cout << "Done" << endl;
	cout << "Time (ms):" << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << endl;
	return 0;
}