//============================================================================
// Name        : keller_miksis_thrust.cu
// Author      : Lambert Plavecz
// Version     :
// Copyright   : no
// Description : Parameter study of the Keller-Miksis equation with odeint RKCK54 and Thrust
//============================================================================


#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <chrono>
#define _USE_MATH_DEFINES
#include <cmath>
#include <utility>

#include <thrust/device_vector.h>
#include <thrust/reduce.h>
#include <thrust/functional.h>

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/thrust/thrust.hpp>


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

typedef double value_type;
typedef thrust::device_vector< value_type > state_type;
typedef thrust::device_vector< int > int_vector;

const int num = 128; 

string file_name = "kellermiksis_thrust_output.txt";

struct keller_miksis
{
	size_t m_N;
    state_type &m_f;
	state_type C;
	
	struct param_fill_functor
    {
		
        template< class T >
        __host__ __device__
        void operator()( T t )const
        {	
            value_type omega = thrust::get< 1 >( t );
			double twr = 2*M_PI/(R_E*omega);
			double P_A1 = 1.5e5;
			double omega2 = 0.0;
			double P_A2 = 0.0;
            thrust::get< 0 >( thrust::get< 0 >( t )) = (P_inf - p_v + 2*sigma/R_E)/ro_L*twr*twr;
			thrust::get< 1 >( thrust::get< 0 >( t )) = (1-3*gam)/(ro_L*c_L)*(P_inf - p_v + 2*sigma/R_E)*twr; 
			thrust::get< 2 >( thrust::get< 0 >( t )) = (P_inf - p_v)/ro_L * twr*twr;
			thrust::get< 3 >( thrust::get< 0 >( t )) = 2*sigma/(ro_L*R_E) *twr*twr;
			thrust::get< 4 >( thrust::get< 0 >( t )) = 4*mu_L/(ro_L*R_E*R_E) * 2*M_PI/omega;
			thrust::get< 5 >( thrust::get< 0 >( t )) = P_A1/ro_L * twr*twr;
			thrust::get< 6 >( thrust::get< 0 >( t )) = P_A2/ro_L *twr*twr;
			thrust::get< 0 >( thrust::get< 2 >( t )) = R_E * omega*P_A1/(ro_L*c_L) * twr*twr;
			thrust::get< 1 >( thrust::get< 2 >( t )) = R_E * omega*P_A2/(ro_L*c_L) * twr*twr;
			thrust::get< 2 >( thrust::get< 2 >( t )) = R_E*omega/(2*M_PI*c_L);
			thrust::get< 3 >( thrust::get< 2 >( t )) = 3*gam;
			thrust::get< 4 >( thrust::get< 2 >( t )) = omega2/omega;
			thrust::get< 5 >( thrust::get< 2 >( t )) =  theta;
        }
    };

    keller_miksis( size_t N , state_type &f): m_N( N ) , m_f( f )
	{ 
		C = state_type(m_N*13); //create param vector
		//fill param vector
		thrust::for_each(
                thrust::make_zip_iterator( thrust::make_tuple(
                        thrust::make_zip_iterator( thrust::make_tuple(
							C.begin(),
							C.begin() + m_N ,
							C.begin() + 2 * m_N,
							C.begin() + 3 * m_N ,
							C.begin() + 4 * m_N , 
							C.begin() + 5 * m_N , 
							C.begin() + 6 * m_N 
						)),
                        m_f.begin() ,
                        thrust::make_zip_iterator( thrust::make_tuple(
							C.begin() + 7 * m_N,
							C.begin() + 8 * m_N ,
							C.begin() + 9 * m_N , 
							C.begin() + 10 * m_N , 
							C.begin() + 11 * m_N ,
							C.begin() + 12 * m_N 
						))
				) ) ,
                thrust::make_zip_iterator( thrust::make_tuple(
                        thrust::make_zip_iterator( thrust::make_tuple(
							C.begin() + m_N ,
							C.begin() + 2 * m_N,
							C.begin() + 3 * m_N ,
							C.begin() + 4 * m_N , 
							C.begin() + 5 * m_N , 
							C.begin() + 6 * m_N ,
							C.begin() + 7 * m_N 
						)),
                        m_f.end() ,
                        thrust::make_zip_iterator( thrust::make_tuple(
							C.begin() + 8 * m_N ,
							C.begin() + 9 * m_N , 
							C.begin() + 10 * m_N , 
							C.begin() + 11 * m_N ,
							C.begin() + 12 * m_N ,
							C.end()
						))
				) ) ,
                param_fill_functor() );
		//f is no longer needed, free it
		m_f.clear();
		m_f.shrink_to_fit();
	}
	
	struct impact_functor
    {
		
		double m_time;
		impact_functor(double time): m_time(time){}
		
        template< class T >
        __host__ __device__
        void operator()( T t )const
        {
            value_type q = thrust::get< 4 >( t );
            value_type x = thrust::get< 0 >( t );
            value_type y = thrust::get< 1 >( t );
			value_type C0 = thrust::get< 0 >( thrust::get< 5 >( t ));
			value_type C1 = thrust::get< 1 >( thrust::get< 5 >( t ));
			value_type C2 = thrust::get< 2 >( thrust::get< 5 >( t ));
			value_type C3 = thrust::get< 3 >( thrust::get< 5 >( t ));
			value_type C4 = thrust::get< 4 >( thrust::get< 5 >( t ));
			value_type C5 = thrust::get< 5 >( thrust::get< 5 >( t ));
			value_type C6 = thrust::get< 6 >( thrust::get< 5 >( t )); 
			value_type C7 = thrust::get< 0 >( thrust::get< 6 >( t ));
			value_type C8 = thrust::get< 1 >( thrust::get< 6 >( t ));
			value_type C9 = thrust::get< 2 >( thrust::get< 6 >( t ));
			value_type C10 = thrust::get< 3 >( thrust::get< 6 >( t ));
			value_type C11 = thrust::get< 4 >( thrust::get< 6 >( t ));
			value_type C12 = thrust::get< 5 >( thrust::get< 6 >( t ));
			
			double rx0 = 1.0/x;
			double N = (C0+C1*y)*pow(rx0,C10) - C2*(1.0+C9*y) -C3*rx0 -C4*y*rx0 - (1.0 - C9*y/3.0)*1.5*y*y
						-(C5*sin(2.0*M_PI*m_time) + C6*sin(2.0*M_PI*C11*m_time + C12)) * (1.0+C9*y)
						-x*(C7*cos(2.0*M_PI*m_time) + C8*cos(2.0*M_PI*C11*m_time+C12));

			double D = x - C9*x*y + C4*C9;
			
            thrust::get< 2 >( t ) = y;
            thrust::get< 3 >( t ) = N/D;
        }
    };

    template< class State , class Deriv >
    void operator()(  const State &x , Deriv &dxdt , value_type t ) const
    {
        thrust::for_each(
                thrust::make_zip_iterator( thrust::make_tuple(
                        boost::begin( x ) ,
                        boost::begin( x ) + m_N ,
						boost::begin( dxdt ) ,
                        boost::begin( dxdt ) + m_N ,
                        m_f.begin() ,
						thrust::make_zip_iterator( thrust::make_tuple(
							C.begin(),
							C.begin() + m_N ,
							C.begin() + 2 * m_N,
							C.begin() + 3 * m_N ,
							C.begin() + 4 * m_N , 
							C.begin() + 5 * m_N , 
							C.begin() + 6 * m_N 
						) ),
						thrust::make_zip_iterator( thrust::make_tuple(
							C.begin() + 7 * m_N,
							C.begin() + 8 * m_N ,
							C.begin() + 9 * m_N , 
							C.begin() + 10 * m_N , 
							C.begin() + 11 * m_N ,
							C.begin() + 12 * m_N 
						) )
				) ),
                thrust::make_zip_iterator( thrust::make_tuple(
                        boost::begin( x ) + m_N ,
                        boost::begin( x ) + 2 * m_N ,
						boost::begin( dxdt ) + m_N ,
                        boost::begin( dxdt ) + 2 * m_N ,
                        m_f.end() ,
						thrust::make_zip_iterator( thrust::make_tuple(
							C.begin() + m_N ,
							C.begin() + 2 * m_N,
							C.begin() + 3 * m_N ,
							C.begin() + 4 * m_N , 
							C.begin() + 5 * m_N , 
							C.begin() + 6 * m_N ,
							C.begin() + 7 * m_N 
						) ),
						thrust::make_zip_iterator( thrust::make_tuple(
							C.begin() + 8 * m_N ,
							C.begin() + 9 * m_N , 
							C.begin() + 10 * m_N , 
							C.begin() + 11 * m_N ,
							C.begin() + 12 * m_N ,
							C.end()
						) )
                ) ) ,
                impact_functor(t) );
    }
};

class observer
{

public:

	struct observer_functor
    {
		
		template< class T >
        __host__ __device__
        void operator()( T t )
        {
			// unpack the parameter we want to vary and the Lorenz variables
            int count = thrust::get< 3 >( t );
			if(count > 2048+64) return;
            value_type x = thrust::get< 0 >( t );
            value_type y = thrust::get< 1 >( t );
			value_type y_prev = thrust::get< 2 >( t );
			value_type extr_prev = thrust::get< 4 >( t );
			value_type min = thrust::get< 5 >( t );
			value_type max = thrust::get< 6 >( t );
			if(y*y_prev < 0.0){//extremum
				if(count > 2047 && count < 2048 + 64){
					if(x > max) thrust::get< 6 >( t ) = x; //max
					if(x < min) thrust::get< 5 >( t ) = x; //min
				}
				count++;
				
				if(fabs(extr_prev - x) < 1.0e-9 && fabs(x) > 1.0e-9 && y_prev!=0.0){ //convergence detection
					if(x > max) thrust::get< 6 >( t ) = x; //max
					if(x < min) thrust::get< 5 >( t ) = x; //min
					count = 5000;
				}
				thrust::get< 4 >( t ) = x;
				thrust::get< 3 >( t ) = count;
			}
			thrust::get< 2 >( t ) = y;
        }
    };
	
    observer(size_t N, state_type &min, state_type &max, int_vector &count): m_N( N ), m_max(max), m_min(min), m_count(count)
	{
    	y_prev = state_type( N );
		thrust::fill(y_prev.begin(), y_prev.end(), 0.0);
		thrust::fill(m_count.begin(), m_count.end(), 0);
		extr_prev = state_type( N );
		thrust::fill(extr_prev.begin(), extr_prev.end(), 100.0); //arbitrary large number
    }

	template< class State >
    void operator()( State &x, double t )
	//(const state_type &x , value_type t )
	{
		thrust::for_each(
                thrust::make_zip_iterator( thrust::make_tuple(
                        x.begin() ,
                        x.begin() + m_N ,
						y_prev.begin(),
						m_count.begin() ,
						extr_prev.begin(),
						m_min.begin(),
						m_max.begin() ) ),
                thrust::make_zip_iterator( thrust::make_tuple(
                        x.begin() + m_N ,
                        x.end(),
						y_prev.end() ,
						m_count.end() ,
						extr_prev.end(),
						m_min.end() ,
						m_max.end() ) ) ,
                obs_fun );
		
		auto min_pos = thrust::min_element(m_count.begin(),m_count.end())-m_count.begin();
		if(m_count[min_pos] >= 2048+64) throw "Ended";
    }
	
private:
	state_type extr_prev;
	state_type y_prev;
	state_type &m_min;
	state_type &m_max;
	size_t m_N;
	int_vector &m_count;
	observer_functor obs_fun;
};
int nums[18] = {256, 768, 1536, 3072, 3840, 5120, 7680, 15360, 30720, 46080, 61440, 76800, 92160, 122880, 184320, 307200, 768000, 4147200};

int main() {
	cout << "Keller-Miksis Thrust started" << endl;

	typedef runge_kutta_cash_karp54< state_type , value_type , state_type , value_type > stepper_type;

	//for(int jj=0; jj < 18;jj++){ //parameter loop
	
	//num = nums[jj];
	cout << num << endl;
	auto t1 = chrono::high_resolution_clock::now();

	thrust::host_vector< value_type > f_host(num);
	double B = 20.0;
	double E = 1000.0;
	double EpB = E/B;
	double invnum = 1.0/(num-1);
	for( size_t i=0 ; i<num ; i++){
		f_host[i] = B*pow(EpB, i*invnum)*2000.0*M_PI;
	}

	state_type f = f_host;

	state_type x( 2 * num );

	// initialize x
	thrust::fill( x.begin() , x.begin() + num , 1.0 );
	// initialize y
	thrust::fill( x.begin() + num, x.end() , 0.0 );
	
	state_type min(num);
	state_type max(num);
	int_vector count(num);
	thrust::fill(count.begin(), count.end(), 0);
	thrust::fill( min.begin(), min.end(), 10000.0); //arbitrary large number
	thrust::fill( max.begin(), max.end(), 0.0);     //arbitrary small number

	keller_miksis km( num , f );
	observer obs(num, min, max, count);
	
	auto stepper = make_controlled(1.0e-10, 1.0e-10, stepper_type());
	
	try{
		integrate_adaptive(stepper, km, x, 0.0, 1.0e23, 0.01, obs);
	}catch(...){
		//cout << "Enough" << endl;
	}
	thrust::host_vector<value_type> min_host(num);
	min_host = min;
	thrust::host_vector<value_type> max_host(num);
	max_host = max;
	
	ofstream ofs(file_name);
	if(!ofs.is_open())exit(-1);
	ofs << setprecision(17);
	for(int u = 0;u < num;u++){
		ofs << 1.5 << " " << f_host[u]/(2000.0*M_PI) << " " << 0.0 << " " << 0.0 << " " << theta << " " << R_E 
			<< " " << min_host[u] << " " << max_host[u] << "\n";
	}

	ofs.flush();
	ofs.close();

	auto t2 = chrono::high_resolution_clock::now();
	cout << "Done" << endl;
	cout << "Time (ms):" << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << endl;
	
	//} //end of parameter loop
	return 0;
}
