//============================================================================
// Name        : valve_vcl.cpp
// Author      : Lambert Plavecz
// Version     : 
// Copyright   : no
// Description : Parameter study of a pressure relief valve system with odeint rkck54
//					includes basic rootfinding algorithm
//============================================================================


#include <iostream>
#include <fstream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <string>
#include <chrono>
#include <boost/numeric/odeint.hpp>
#include "vcl/version2/vectorclass.h"
#include "vcl/version2/vectormath_trig.h"
#include "vcl/version2/vectormath_exp.h"

using namespace boost::numeric::odeint;

const double kappa = 1.25;
const double beta = 20.0;
const double delta = 10.0;
const double r = 0.8;

const double root_tolerance = 1.0e-6;

int num = 256;
const int vectorSize = 4; //number of doubles in vector_type
const int arraySize = 3; //number of equations in ode system

std::string file_name = "impact_dyn_vcl_output.txt";

typedef Vec4d vector_type;
typedef std::array< vector_type , arraySize> state_type;

std::vector<double> extrema(64 * num);


class custom_error_checker
{
public:

    typedef array_algebra algebra_type;
    typedef default_operations operations_type;

    custom_error_checker(
            double eps_abs = static_cast< double >( 1.0e-10 ) ,
            double eps_rel = static_cast< double >( 1.0e-10 ) ,
            double a_x = static_cast< double >( 1.0 ) ,
            double a_dxdt = static_cast< double >( 1.0 ) )
    : m_eps_abs( eps_abs ) , m_eps_rel( eps_rel ) ,
      m_a_x( a_x ) , m_a_dxdt( a_dxdt )
    { }


    template< class State , class Deriv , class Err , class Time >
    double error( const State &x_old ,
                      const Deriv &dxdt_old ,
                      Err &x_err , Time dt ) const
    {
        return error( algebra_type() , x_old , dxdt_old , x_err , dt );
    }

    template< class State , class Deriv , class Err , class Time >
    double error( algebra_type &algebra ,
                      const State &x_old ,
                      const Deriv &dxdt_old ,
                      Err &x_err , Time dt ) const
    {
        // this overwrites x_err 
        algebra.for_each3( x_err , x_old , dxdt_old ,
            typename operations_type::template rel_error< double >(
                m_eps_abs , m_eps_rel , m_a_x ,
                m_a_dxdt * get_unit_value( dt ) ) );

        state_type temp = x_err;
		//take absolute value of every element
		for(int j=0; j < arraySize;j++){
			temp[j] = abs(x_err[j]); //take absolute value of vcl vector elements
		}
		Vec4d v = 0.0;
		//find vcl vector of maximal elements
		for(int j=0; j < arraySize;j++){
			v = max(v,temp[j]);
		}
		//find maximum of elments in maximal vector
		double m = v[0];
		for(int j=1; j < vectorSize;j++){
			if(m < v[j]) m = v[j];
		}
		return m; //return greatest element
    }

private:

    double m_eps_abs;
    double m_eps_rel;
    double m_a_x;
    double m_a_dxdt;
};

class valve{
	vector_type &m_q;

public:
	int system_number = 0;

	valve(vector_type &q): m_q(q)
	{}
	void operator() (const state_type &x, state_type &dxdt, const double t){
		//x[0] displacement, x[1] velocity, x[2] pressure in chamber
		dxdt[0] = x[1];
    	dxdt[1] = -kappa*x[1] - (x[0] + delta) + x[2];
    	dxdt[2] = beta * (m_q - x[0] * sqrt(x[2]));
	}
	void operator() (const std::array<double, arraySize> &x, std::array<double, arraySize> &dxdt, double t){ 
		dxdt[0] = x[1];
    	dxdt[1] = -kappa*x[1] - (x[0] + delta) + x[2];
    	dxdt[2] = beta * (m_q[system_number] - x[0]* std::sqrt(x[2]));
	}
};

using namespace std;

int min_int_array(array<int, vectorSize> t){
	int m = t[0];
	for(int i=1; i < vectorSize;i++){
		if(m > t[i]) m = t[i];
	}
	return m;
}

typedef runge_kutta_cash_karp54< array<double, arraySize>, double , array<double, arraySize>, double > observer_stepper;

class observer
{
    state_type x_prev;
	double t_prev;
    array<int, vectorSize> count;
	array<double, vectorSize> extr_prev;
	double * tmp_array;
    int no; //id number of ode
	valve & m_system; //the ode system
	controlled_runge_kutta<observer_stepper> m_stepper = make_controlled( 1e-10, 1e-10, observer_stepper());

public:
    observer(int row, valve & equations): no(row), m_system(equations) 
	{ 	
		tmp_array = new double[vectorSize*arraySize];
		for(int j=0; j < vectorSize;j++){
			extr_prev[j] = 200.0;
			count[j] = 0;
		}
    }
	~observer()
	{
		delete[] tmp_array;
	}

    void operator()(state_type &x , double t )
    {
    	if(t==0.0){
    		x_prev = x;
			t_prev = t;
    		return;
    	}
		
		for(int i=0; i < arraySize;i++)
			for(int j=0; j < vectorSize; j++)
				tmp_array[i*vectorSize+j] = x[i][j]; //fill tmp_array with the current x values
			
		for(int j=0; j < vectorSize;j++){
			if(fabs(x[0][j]) < root_tolerance){
				if(count[j] > 2047 && count[j] < 2048 + 64){ //saving
					extrema[64*(no+j)+count[j]-2048] = x[0][j];
				}else{ //extremum search ended in this system
					if(min_int_array(count) >= 2048 + 64) throw 0;
				}
				count[j]++;
				if(count[j] < 2048+64 && fabs(extr_prev[j] - x[0][j]) < 1.0e-9 && fabs(x[0][j]) > 1.0e-9){ //check whether solution converges
					count[j] = 5000;
					//if true end integration
					if(min_int_array(count) >= 2048 + 64) throw 0;
				}
				extr_prev[j] = x[0][j];
				tmp_array[1*vectorSize+j] = -r * x[1][j]; //impact
			}
			else if(fabs(x[1][j]) < root_tolerance){
				if(count[j] > 2047 && count[j] < 2048 + 64){ //saving
					extrema[64*(no+j)+count[j]-2048] = x[0][j];
				}else{ //extremum search ended in this system
					if(min_int_array(count) >= 2048 + 64) throw 0;
				}
				count[j]++;
				if(count[j] < 2048+64 && fabs(extr_prev[j] - x[0][j]) < 1.0e-9 && fabs(x[0][j]) > 1.0e-9){ //check whether solution converges
					count[j] = 5000;
					if(min_int_array(count) >= 2048 + 64) throw 0;
				}
				extr_prev[j] = x[0][j];
			}
			else if(x[0][j]*x_prev[0][j] < 0.0 && x[1][j] < 0.0){//extremum
				if(count[j] < 2048 + 64){ //saving
					
					double t_start = t_prev;
					double t_start_tmp = t_start;
					double t_end = t;
					array<double, arraySize> x_start;
					x_start[0] = x_prev[0][j];
					x_start[1] = x_prev[1][j];
					x_start[2] = x_prev[2][j];
					array<double, arraySize> x_tmp = x_start;
				
					double dt = (t_end-t_start)/2;
					double dt_temp = dt;
					m_system.system_number = j;
					m_stepper.try_step(boost::ref(m_system), x_tmp, t_start_tmp, dt_temp);
					while(fabs(x_tmp[0]) > root_tolerance && dt!=0.0){
						if(x_start[0]*x_tmp[0] < 0.0){
							dt = dt/2;
							dt_temp = dt;
							t_start_tmp = t_start;
							x_tmp = x_start;
							m_stepper.try_step(boost::ref(m_system), x_tmp, t_start_tmp, dt_temp);
						}else{
							t_start += dt;
							dt = dt/2;
							dt_temp = dt;
							t_start_tmp = t_start;
							x_start = x_tmp;
							m_stepper.try_step(boost::ref(m_system), x_tmp, t_start_tmp, dt_temp);
						}
					}
					//x_tmp is the array of the event x coordinates
					if(count[j] > 2047) extrema[64*(no+j)+count[j]-2048] = x_tmp[0];
										
					x_tmp[1] = -r * x_tmp[1]; //impact
					
					dt_temp =  t-t_start_tmp;
					m_stepper.try_step(boost::ref(m_system), x_tmp, t_start_tmp, dt_temp); //one step to t
					//x_tmp: corrected values at t time, insert them into tmp_array
					for(int i=0; i < arraySize; i++)
						tmp_array[i*vectorSize+j] = x_tmp[i];
					
				}else{ //extremum search ended
					if(min_int_array(count) > 2048 + 63) throw 0;
				}
				count[j]++;
				if(count[j] < 2048+64 && fabs(extr_prev[j] - x[0][j]) < 1.0e-9 && fabs(x[0][j]) > 1.0e-9){ //check whether solution converges
					//if true end integration
					count[j] = 5000;
					if(min_int_array(count) > 2048 + 63) throw 0;
				}
				extr_prev[j] = x[0][j];
			}
			else if(x[1][j]*x_prev[1][j] < 0.0 && fabs(x_prev[1][j]) > root_tolerance){//extremum
				if(count[j] > 2047 && count[j] < 2048 + 64){ //saving
					
					double t_start = t_prev;
					double t_start_tmp = t_start;
					double t_end = t;
					array<double, arraySize> x_start;
					x_start[0] = x_prev[0][j];
					x_start[1] = x_prev[1][j];
					x_start[2] = x_prev[2][j];
					array<double, arraySize> x_tmp = x_start;
				
					double dt = (t_end-t_start)/2;
					double dt_temp = dt;
					m_system.system_number = j;
					m_stepper.try_step(boost::ref(m_system), x_tmp, t_start_tmp, dt_temp);
					while(fabs(x_tmp[0]) > root_tolerance && dt!=0.0){
						if(x_start[0]*x_tmp[0] < 0.0){
							dt = dt/2;
							dt_temp = dt;
							t_start_tmp = t_start;
							x_tmp = x_start;
							m_stepper.try_step(boost::ref(m_system), x_tmp, t_start_tmp, dt_temp);
						}else{
							t_start += dt;
							dt = dt/2;
							dt_temp = dt;
							t_start_tmp = t_start;
							x_start = x_tmp;
							m_stepper.try_step(boost::ref(m_system), x_tmp, t_start_tmp, dt_temp);
						}
						
					}
					//x_tmp is the array of the event x coordinates
					extrema[64*(no+j)+count[j]-2048] = x[0][j];
					
					dt_temp =  t-t_start_tmp;
					m_stepper.try_step(boost::ref(m_system), x_tmp, t_start_tmp, dt_temp); //one step to t
					//x_tmp: corrected values at t time, insert them into tmp_array
					for(int i=0; i < arraySize; i++)
						tmp_array[i*vectorSize+j] = x_tmp[i];
				}else{ //extremum search ended
					if(min_int_array(count) > 2048 + 63) throw 0;
				}
				count[j]++;
				if(count[j] < 2048+64 && fabs(extr_prev[j] - x[0][j]) < 1.0e-9 && fabs(x[0][j]) > 1.0e-9){ //check whether solution converges
					//if true save the current value and end integration
					count[j] = 5000;
					if(min_int_array(count) > 2048 + 63) throw 0;
				}
				extr_prev[j] = x[0][j];
			}
		}
		for(int i=0; i < arraySize;i++){
			x[i].load(tmp_array + i*vectorSize); //copy the corrected tmp_array back into the state_type
		}
		x_prev = x;
		t_prev = t;
    }
};


int nums[11] = {256, 768, 1536, 3072, 3840, 5120, 	7680,    15360,   30720, 46080, 61440}; // 76800, 92160, 122880, 184320, 307200, 768000, 4147200};
//			 11029 27909 52764 100145 125182 
int main() {
	cout << "Begin" << endl;
	
	typedef runge_kutta_cash_karp54< state_type , double , state_type , double > stepper_type;
	typedef controlled_runge_kutta< stepper_type, custom_error_checker > controlled_stepper;
	
	controlled_stepper stepper(custom_error_checker(1.0e-10, 1.0e-10));
	state_type x;
	
	double * temp_arr = new double[vectorSize];
	
for(int jj=0; jj < 11;jj++){ //parameter loop
	
	num = nums[jj];
	cout << num << endl;
	
	extrema.resize(64 * num);

	auto t1 = chrono::high_resolution_clock::now();

	double invnum = 1.0/(num-1);
	vector_type q = 0.0;
	for(int u = 0;u < num;u+=vectorSize){ //f1
		//cout << u << endl;		
		for(int i=0;i < vectorSize;i++){ 
			temp_arr[i] = 0.2 + 9.8*(u+i)*invnum;
		}			
		
		x[0] = 0.2;
		x[1] = 0.0;
		x[2] = 0.0;
		
		q.load(temp_arr);
		valve equation(q);
		observer obs(u, equation);

		try{
			integrate_adaptive(stepper, boost::ref(equation), x, 0.0, 1e23, 0.01, boost::ref(obs));
		}catch(...){
			//cout << "Enough" << endl;
		}
	}

	ofstream ofs(file_name);
	if(!ofs.is_open())exit(-1);
	ofs.precision(17);
	//ofs.flags(ios::scientific);
	for(int u = 0;u < num;u++){
		ofs << 0.2 + 9.8*u*invnum;
		for(int i=0;i < 64;i++){
			ofs << " " << extrema[64*u + i];
		}
		ofs << "\n";
	}

	ofs.flush();
	ofs.close();

	auto t2 = chrono::high_resolution_clock::now();
	cout << "Done" << endl;
	cout << "Time (ms):" << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << endl;
	
}
	delete[] temp_arr;
	
	return 0;
}
