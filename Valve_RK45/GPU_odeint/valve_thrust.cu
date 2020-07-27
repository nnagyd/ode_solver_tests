//============================================================================
// Name        : valve_thrust.cu
// Author      : Lambert Plavecz
// Version     :
// Copyright   : no
// Description : Simulation of a pressure relief valve with odeint CUDA Thrust RKCK54
//============================================================================


#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <chrono>
#include <cmath>
#include <utility>

#include <thrust/device_vector.h>
#include <thrust/reduce.h>
#include <thrust/functional.h>

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/thrust/thrust.hpp>


using namespace std;
using namespace boost::numeric::odeint;

typedef double value_type;
typedef thrust::device_vector< value_type > state_type;
typedef thrust::device_vector< int > int_vector;

int num = 12; //12 - 5952528 ms
value_type inv_num = 1.0/(num-1);

string file_name = "impact_dyn_output_thrust.txt";

const value_type kappa = 1.25;
const value_type beta = 20.0;
const value_type delta = 10.0;
const value_type r = 0.8;


struct impact_dynamics
{
    struct impact_functor
    {
        template< class T >
        __host__ __device__
        void operator()( T t )const
        {												
            value_type q = thrust::get< 3 >( t );
            value_type x = thrust::get< 0 >( t );
            value_type y = thrust::get< 1 >( t );
            value_type z = thrust::get< 2 >( t );
            thrust::get< 4 >( t ) = y;
            thrust::get< 5 >( t ) = -kappa * y - (x + delta) + z;
            thrust::get< 6 >( t ) = beta * (q - x * sqrt(z));

        }
    };

    impact_dynamics( size_t N , const state_type &q)
    : m_N( N ) , m_q( q ) { }

    template< class State , class Deriv >
    void operator()(  const State &x , Deriv &dxdt , value_type t ) const
    {
        thrust::for_each(
                thrust::make_zip_iterator( thrust::make_tuple(
                        x.begin() ,
                        x.begin() + m_N ,
                        x.begin() + 2 * m_N ,
                        m_q.begin() ,
                        dxdt.begin() ,
                        dxdt.begin() + m_N ,
                        dxdt.begin() + 2 * m_N  ) ) ,
                thrust::make_zip_iterator( thrust::make_tuple(
                        x.begin() + m_N ,
                        x.begin() + 2 * m_N ,
                        x.end() ,
                        m_q.end() ,
                        dxdt.begin() + m_N ,
                        dxdt.begin() + 2 * m_N ,
                        dxdt.end()  ) ) ,
                impact_functor() );
    }
    size_t m_N;
    const state_type &m_q;
};

class impact_observer
{
public:

	struct observer_functor
    {
		
		template< class T >
        __host__ __device__
        void operator()( T t )
        {
			int count = thrust::get< 3 >( t );
			if(count > 2048+63)return;
			value_type extr_prev = thrust::get< 4 >( t );
            value_type x = thrust::get< 0 >( t );
            value_type y = thrust::get< 1 >( t );
			value_type y_prev = thrust::get< 2 >( t );
			int extr_place = thrust::get< 5 >( t );
			if(x <= 0 && y < 0){
				y = -r * y;
				thrust::get< 1 >( t ) = y;
			}
			if(y*y_prev < 0){//extremum
				if(2047 < count){
					thrust::get< 6 >( t ) = x;
					thrust::get< 5 >( t ) = extr_place+1;
				}
				count++;
				if(fabs(extr_prev - x) < 1.0e-6 && fabs(x) > 1.0e-6 && y_prev!=0.0){ //convergence detection
					thrust::get< 6 >( t ) = x;
					count = 5000;
				}
				thrust::get< 4 >( t ) = x;
				thrust::get< 3 >( t ) = count;
			}
			thrust::get< 2 >( t ) = y;

        }
		observer_functor(){}
    };
	
    impact_observer(size_t N, state_type &extrema, int_vector &extr_places)
		: m_N( N ), m_extr_places(extr_places), m_extrema(extrema)
	{
		y_prev = state_type(N);
		thrust::fill(y_prev.begin(), y_prev.end(), 0.0);
		m_count = int_vector(N);
		thrust::fill(m_count.begin(), m_count.end(), 0);
		extr_prev = state_type(N);
		thrust::fill(extr_prev.begin(), extr_prev.end(), 10000.0); //arbitrary large number
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
						m_extr_places.begin(),
						thrust::make_permutation_iterator(m_extrema.begin(),m_extr_places.begin())   ) ),
                thrust::make_zip_iterator( thrust::make_tuple(
                        x.begin() + m_N ,
                        x.begin() + 2 * m_N ,
						y_prev.end() ,
						m_count.end() ,
						extr_prev.end(),
						m_extr_places.end(),
						thrust::make_permutation_iterator(m_extrema.begin(),m_extr_places.end())  ) ) ,
                obs_func );
		
		auto min_pos = thrust::min_element(m_count.begin(),m_count.end())-m_count.begin();
		if(min_count < m_count[min_pos]){//cout << m_count[min_pos] << endl;
			min_count = m_count[min_pos];
			if(min_count%200==0) cout << min_count << endl;
			//if( min_count == 1366) for(int z=0;z < m_N; z++) cout << m_count[z] << endl;
		}
		if(m_count[min_pos] > 2048+63) throw 1;
    }
private:	
	state_type y_prev;
	state_type extr_prev;
	state_type &m_extrema; //pointer to device memory
	size_t m_N;
	int min_count = 0;
	int_vector m_count;
	int_vector &m_extr_places;
	observer_functor obs_func;
};
int nums[18] = {256, 768, 1536, 3072, 3840, 5120, 7680, 15360, 30720, 46080, 61440, 76800, 92160, 122880, 184320, 307200, 768000, 4147200};
//343767 370154 389796 410224 398229 250804 453742 0 700243

int main() {
	cout << "Impact dynamics started" << endl;

	//TODO try dopri5
	typedef runge_kutta_cash_karp54< state_type , value_type , state_type , value_type > stepper_type;
	//typedef controlled_runge_kutta< stepper_type, custom_error_checker > controlled_stepper;
	
	/*for(int jj=8; jj < 18; jj++){
		num = nums[jj];
		inv_num = 1.0/num;
		*/
		
	cout << num << endl;
	auto t1 = chrono::high_resolution_clock::now();

	thrust::host_vector< value_type > q_host(num);
	for( int i=0 ; i<num ; i++)
		q_host[i] = 0.2 + value_type(i) * inv_num * 9.8 ;
	
	thrust::host_vector<int> extremum_places_host(num);
	for( int i=0 ; i<num ; i++)
		extremum_places_host[i] = i*64;

	state_type q = q_host;

	state_type x( 3 * num );

	// initialize x
	thrust::fill( x.begin() , x.begin() + num , 0.2 );
	// initialize y,z
	thrust::fill( x.begin() + num, x.end() , 0.0 );
	

	state_type extrema(64 * num);
	thrust::fill( extrema.begin(), extrema.end(), -1.0);
		
	int_vector extremum_places = extremum_places_host;

	impact_dynamics impact( num , q );
	impact_observer obs(num, extrema, extremum_places);
	
	auto stepper = make_controlled( 1.0e-10 , 1.0e-10 , stepper_type());
	
	try{
		integrate_adaptive( stepper , impact, x, 0.0, 1.0e23, 0.01, obs);
	}catch(...){
		//cout << "Enough" << endl;
	}
	thrust::host_vector<value_type> extrema_host(64 * num);
	extrema_host = extrema; //copy back from GPU
	
	ofstream ofs(file_name);
	if(!ofs.is_open())exit(-1);
	ofs.precision(17);
	for(int u = 0;u < num;u++){
		ofs << (0.2 + u * 9.8 * inv_num);
		for(int i=0; i < 64; i++) ofs << " " << extrema_host[u*64+i];
		ofs	<< "\n";
	}

	ofs.flush();
	ofs.close();

	auto t2 = chrono::high_resolution_clock::now();
	cout << "Done" << endl;
	cout << "Time (ms):" << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << endl;
	
	//}

	return 0;
}
