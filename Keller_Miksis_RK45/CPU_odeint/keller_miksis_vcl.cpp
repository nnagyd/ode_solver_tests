//============================================================================
// Name        : keller_miksis_vcl.cpp
// Author      : Lambert Plavecz
// Version     :
// Copyright   : no
// Description : Parameter study of the Keller-Miksis equation with odeint rkck54
//					using explicit vectorisation via Vectorclass Library
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


int num = 0;
const int vectorSize = 4;
const int arraySize = 2;

std::string file_name = "kellermiksis_vcl_noroot_output.txt";

typedef Vec4d value_type;
typedef std::array< value_type , arraySize> state_type;

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


using namespace std;

class km{
	std::array<value_type, 13> C;

public:
	km(value_type P_A1, value_type P_A2, value_type omega, value_type omega2){
		value_type twr = 2*M_PI/(R_E*omega);
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
		value_type rx0 = 1.0/x[0];
		value_type N = (C[0]+C[1]*x[1])*pow(rx0,C[10]) - C[2]*(1.0+C[9]*x[1]) -C[3]*rx0 -C[4]*x[1]*rx0 - (1.0 - C[9]*x[1]/3.0)*1.5*x[1]*x[1]
						-(C[5]*sin(2.0*M_PI*t) + C[6]*sin(2.0*M_PI*C[11]*t + C[12])) * (1.0+C[9]*x[1])
						-x[0]*(C[7]*cos(2.0*M_PI*t) + C[8]*cos(2.0*M_PI*C[11]*t+C[12]));

		value_type D = x[0] - C[9]*x[0]*x[1] + C[4]*C[9];

		dxdt[0] = x[1];
		dxdt[1] = N/D;
	}
};

class observer
{
    int no; //id number of ode

public:
    observer(int row): no(row) { 
    	for(int j=0; j < vectorSize;j++){
			extrema[no+j] = 0.0; //initial max
		}
    }

    void operator()(const state_type &x , double t )
    {
		for(int j=0; j < vectorSize;j++){
			if(x[0][j] > extrema[no+j]) extrema[no+j] = x[0][j]; //max saving
		}

    }
};


int nums[11] = {256, 768, 1536, 3072, 3840, 5120, 	7680,    15360,   30720, 46080, 61440}; // 76800, 92160, 122880, 184320, 307200, 768000, 4147200};
//milliseconds 41633 108161 199804 371821 457157 592988 857944 1633295 3141603
int main() {
	cout << "Begin" << endl;

	typedef runge_kutta_cash_karp54< state_type , double , state_type , double > stepper_type;
	typedef controlled_runge_kutta< stepper_type, custom_error_checker > controlled_stepper ;
	
	controlled_stepper stepper(custom_error_checker(1.0e-10, 1.0e-10));
	state_type x;
	
	double * temp_arr = new double[vectorSize];
	
for(int jj=0; jj < 11;jj++){ //parameter loop
	
	num = nums[jj];
	cout << num << endl;
	extrema = std::vector<double>(vectorSize * num);

	auto t1 = chrono::high_resolution_clock::now();

	double B = 20.0;
	double E = 1000.0;
	double EpB = E/B;
	double invnum = 1.0/(num-1);
	value_type f1 = 0.0;
	for(int u = 0;u < num;u+=vectorSize){ //f1
		//if(u%100==0) cout << u << endl;		
		for(int i=0;i < vectorSize;i++){ 
			temp_arr[i] = B*std::pow(EpB, (u+i)*invnum);
		}			

		observer obs(u);
		
		x[0] = 1.0;
		x[1] = 0.0;
		f1.load(temp_arr);
		
		km equation( 1.5*1e5, 0.0, 2000.0*M_PI*f1, 0.0);

		integrate_adaptive(boost::ref(stepper), boost::ref(equation), x, 0.0, 1024.0, 0.01);
		
		integrate_adaptive(boost::ref(stepper), boost::ref(equation), x, 1024.0, 1088.0, 0.01, boost::ref(obs));
	}

	ofstream ofs(file_name);
	if(!ofs.is_open())exit(-1);
	ofs.precision(17);
	//ofs.flags(ios::scientific);
	for(int u = 0;u < num;u++){
		ofs << 1.5 << " " << B*pow(EpB, u*invnum) << " " << 0.0 << " " << 0.0 << " " << theta << " " << R_E
			<< " " << extrema[u] << "\n";
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
