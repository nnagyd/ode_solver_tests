//============================================================================
// Name        : lorenz_RK4.cpp
// Author      : Lambert Plavecz
// Version     :
// Copyright   : no
// Description : lorenz with odeint CUDA Thrust RK4, Based on the example available at: 
// https://github.com/headmyshoulder/odeint-v2/blob/master/examples/thrust/lorenz_parameters.cu
//============================================================================


#include <iostream>
#include <ctime>
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

int num = 0;
value_type inv_num = 1.0/(num-1);
const int step_number = 1000;

const value_type sigma = 10.0;
const value_type b = 8.0/3.0;
const value_type dt = 0.01;

struct lorenz_system
{
    struct lorenz_functor
    {
        template< class T >
        __host__ __device__
        void operator()( T t )const
        {
            // unpack the parameter we want to vary and the Lorenz variables
            value_type k = thrust::get< 3 >( t );
            value_type x = thrust::get< 0 >( t );
            value_type y = thrust::get< 1 >( t );
            value_type z = thrust::get< 2 >( t );
            thrust::get< 4 >( t ) = sigma * ( y - x );
            thrust::get< 5 >( t ) = k * x - y - x * z;
            thrust::get< 6 >( t ) = -b * z + x * y ;

        }
    };

    lorenz_system( size_t N , const state_type &k )
    : m_N( N ) , m_k( k ) { }

    template< class State , class Deriv >
    void operator()(  const State &x , Deriv &dxdt , value_type t ) const
    {
        thrust::for_each(
                thrust::make_zip_iterator( thrust::make_tuple(
                        x.begin(),
                        x.begin() + m_N ,
                        x.begin() + 2 * m_N ,
                        m_k.begin() ,
                        dxdt.begin() ,
                        dxdt.begin() + m_N ,
                        dxdt.begin() + 2 * m_N  ) ) ,
                thrust::make_zip_iterator( thrust::make_tuple(
                        x.begin() + m_N ,
                        x.begin() + 2 * m_N ,
                        x.end(),
                        m_k.end() ,
                        dxdt.begin() + m_N ,
                        dxdt.begin() + 2 * m_N ,
                        dxdt.end() ) ) ,
                lorenz_functor() );
    }
    size_t m_N;
    const state_type &m_k;
};

int nums[13] = {256, 512,  1024,  2048,  4096, 7680, 15360, 46080, 92160, 184320, 307200, 768000, 4147200};
int main() {
	cout << "Lorenz RK4 started" << endl;

	runge_kutta4<state_type> stepper;
	
	for(int j=0; j < 13;j++){
		num = nums[j];
		inv_num = 1.0/(num-1);
		cout << num << endl;
		auto t1 = clock();

		vector< value_type > k_host(num);
		for( size_t i=0 ; i<num ; i++)
			k_host[i] = value_type(i) * inv_num * 21.0 ;

		state_type k = k_host;

		state_type x( 3 * num );

		// initialize x,y,z
		thrust::fill( x.begin() , x.end() , 10.0 );

		lorenz_system lorenz( num , k );
	
		integrate_const(stepper, lorenz, x, 0.0, step_number*dt, dt);


		auto t2 = clock();
		cout << "Time (ms):" << (t2 - t1)*1000/CLOCKS_PER_SEC << endl;
	}

	return 0;
}
