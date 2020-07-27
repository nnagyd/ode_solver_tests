//============================================================================
// Name        : keller_miksis.cpp
// Author      : Lambert Plavecz
// Version     : 
// Copyright   : no
// Description : Parameter study of the Keller-Miksis equation with odeint rkck54
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

int num = 256;
const int vectorSize = 4;
const int arraySize = 2;

std::string file_name = "kellermiksis_vcl_output.txt";

typedef Vec4d vector_type;
typedef std::array< vector_type , arraySize> state_type;

std::vector<double> extrema(2 * num);


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
			//take absolute value of vcl vector elements
			temp[j] = abs(x_err[j]);
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

using namespace std;

class km{
	std::array<vector_type, 13> C;

public:
	int system_number = 0;

	km(vector_type P_A1, vector_type P_A2, vector_type omega, vector_type omega2){
		vector_type twr = 2*M_PI/(R_E*omega);
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
		vector_type rx0 = 1.0/x[0];
		vector_type N = (C[0]+C[1]*x[1])*pow(rx0,C[10]) - C[2]*(1.0+C[9]*x[1]) -C[3]*rx0 -C[4]*x[1]*rx0 - (1.0 - C[9]*x[1]/3.0)*1.5*x[1]*x[1]
						-(C[5]*sin(2.0*M_PI*t) + C[6]*sin(2.0*M_PI*C[11]*t + C[12])) * (1.0+C[9]*x[1])
						-x[0]*(C[7]*cos(2.0*M_PI*t) + C[8]*cos(2.0*M_PI*C[11]*t+C[12]));

		vector_type D = x[0] - C[9]*x[0]*x[1] + C[4]*C[9];

		dxdt[0] = x[1];
		dxdt[1] = N/D;
	}
	void operator() (const array<double, 2> &x, array<double, 2> &dxdt, double t){ //FIXME fix me 
		double rx0 = 1.0/x[0];
		double N = (C[0][system_number]+C[1][system_number]*x[1])*pow(rx0,C[10][system_number]) - C[2][system_number]*(1.0+C[9][system_number]*x[1]) -C[3][system_number]*rx0 -C[4][system_number]*x[1]*rx0 - (1.0 - C[9][system_number]*x[1]/3.0)*1.5*x[1]*x[1]
						-(C[5][system_number]*sin(2.0*M_PI*t) + C[6][system_number]*sin(2.0*M_PI*C[11][system_number]*t + C[12][system_number])) * (1.0+C[9][system_number]*x[1])
						-x[0]*(C[7][system_number]*cos(2.0*M_PI*t) + C[8][system_number]*cos(2.0*M_PI*C[11][system_number]*t+C[12][system_number]));

		double D = x[0] - C[9][system_number]*x[0]*x[1] + C[4][system_number]*C[9][system_number];

		dxdt[0] = x[1];
		dxdt[1] = N/D;
	}
};

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
	km & m_system; //the ode system
	controlled_runge_kutta<observer_stepper> m_stepper = make_controlled( 1e-10, 1e-10, observer_stepper());

public:
    observer(int row, km & equations): no(row), m_system(equations) 
	{ 	
		tmp_array = new double[vectorSize*arraySize];
		for(int j=0; j < vectorSize;j++){
			extrema[2*(no+j)] = 0.0; //initial max
			extrema[2*(no+j)+1] = 10000.0; //initial min
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
		bool event = false;
		for(int j=0; j < vectorSize;j++){ //check whether a root finding is necessary
			if(x[1][j]*x_prev[1][j] < 0 && fabs(x_prev[1][j]) > root_tolerance){
				event = true;
				break;
			}
		}
		if(event){
			for(int i=0; i < arraySize;i++)
				for(int j=0; j < vectorSize; j++)
					tmp_array[i*vectorSize+j] = x[i][j]; //fill tmp_array with the current x values
		}
		for(int j=0; j < vectorSize;j++){
			
			if(fabs(x[1][j]) < root_tolerance){
				if(count[j] > 2047 && count[j] < 2048 + 64){ //saving
					if(x[0][j] > extrema[2*(no+j)]) extrema[2*(no+j)] = x[0][j]; //max
					if(x[0][j] < extrema[2*(no+j)+1]) extrema[2*(no+j)+1] = x[0][j]; //min
				}else{ //extremum search ended in this system
					if(min_int_array(count) >= 2048 + 64) throw 0;
				}
				count[j]++;
				if(count[j] < 2048+64 && fabs(extr_prev[j] - x[0][j]) < 1.0e-9 && fabs(x[0][j]) > 1.0e-9){ //check whether solution converges
					//if true save the current value and end integration
					if(x[0][j] > extrema[2*(no+j)]) extrema[2*(no+j)] = x[0][j]; //max
					if(x[0][j] < extrema[2*(no+j)+1]) extrema[2*(no+j)+1] = x[0][j]; //min
					count[j] = 5000;
					if(min_int_array(count) >= 2048 + 64) throw 0;
				}
				extr_prev[j] = x[0][j];
			}
			else if(x[1][j]*x_prev[1][j] < 0.0 && fabs(x_prev[1][j]) > root_tolerance){//extremum
				if(count[j] > 2047 && count[j] < 2048 + 64){ //saving
					
					double t_start = t_prev;
					double t_start_tmp = t_start;
					double t_end = t;
					array<double, 2> x_start;
					x_start[0] = x_prev[0][j];
					x_start[1] = x_prev[1][j];
					array<double, 2> x_tmp = x_start;
				
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
					if(x_tmp[0] > extrema[2*(no+j)]) extrema[2*(no+j)] = x_tmp[0]; //max
					if(x_tmp[0] < extrema[2*(no+j)+1]) extrema[2*(no+j)+1] = x_tmp[0]; //min
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
					//cout << "x_prev: " << x_prev << ", y_prev:" << y_prev << ", x: " << x[0] << endl;
					//if true save the current value and end integration
					if(x[0][j] > extrema[2*(no+j)]) extrema[2*(no+j)] = x[0][j]; //max
					if(x[0][j] < extrema[2*(no+j)+1]) extrema[2*(no+j)+1] = x[0][j]; //min
					count[j] = 5000;
					if(min_int_array(count) > 2048 + 63) throw 0;
				}
				extr_prev[j] = x[0][j];
			}
		}
		if(event) {
			for(int i=0; i < arraySize;i++){
				x[i].load(tmp_array + i*vectorSize); //copy the corrected tmp_array back into the state_type
			}
		}
		x_prev = x;
		t_prev = t;
    }
};


int nums[11] = {256, 768, 1536, 3072, 3840, 5120, 	7680,    15360,   30720, 46080, 61440}; // 76800, 92160, 122880, 184320, 307200, 768000, 4147200};
//			 52629 146096 278676 532625 658181 883737 1259082
int main() {
	cout << "Begin" << endl;
	
	typedef runge_kutta_cash_karp54< state_type , double , state_type , double > stepper_type;
	typedef controlled_runge_kutta< stepper_type, custom_error_checker > controlled_stepper;
	
	controlled_stepper stepper(custom_error_checker(1.0e-10, 1.0e-10));
	state_type x;
	
	double * temp_arr = new double[vectorSize];
	
for(int jj=4; jj < 11;jj++){ //parameter loop
	
	num = nums[jj];
	cout << num << endl;
	extrema = std::vector<double>(2 * num);

	auto t1 = chrono::high_resolution_clock::now();

	double B = 20.0;
	double E = 1000.0;
	double EpB = E/B;
	double invnum = 1.0/(num-1);
	vector_type f1 = 0.0;
	for(int u = 0;u < num;u+=vectorSize){ //f1
		//if(u%100==0) cout << u << endl;		
		for(int i=0;i < vectorSize;i++){ 
			temp_arr[i] = B*std::pow(EpB, (u+i)*invnum);
		}			
		
		x[0] = 1.0;
		x[1] = 0.0;
		f1.load(temp_arr);
		
		km equation( 1.5*1e5, 0.0, 2000.0*M_PI*f1, 0.0);
		observer obs(u, equation);

		try{
			integrate_adaptive(stepper, boost::ref(equation), x, 0.0, 1.0e23, 0.01, boost::ref(obs));
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
	delete[] temp_arr;
	
	return 0;
}
