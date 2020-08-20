### Integrator settings
* Number of steps = 1000
* Number of parameters = 2^17 = 131 072

### Results
|                          |                     |              |               |                  |            | 
|--------------------------|---------------------|--------------|---------------|------------------|------------| 
| "Code"                   | "simple_vcl_unroll" | "simple_vcl" | "simple_loop" | "transcendental" | "division" | 
| "Unroll"                 | 8.                  | 1.           | 1.            | 2.               | 2.         | 
| "Runtimes"               | 0.1342              | 0.5975       | 2.383         | 2.125            | 1.004      | 
| "Deviation"              | 0.0002              | 0.001        | 0.0047        | 0.1645           | 0.0018     | 
| "Clock Cycles"           | 4.96e8              | 2.209e9      | 8.811e9       | 7.643e9          | 3.711e9    | 
| "Number of Instructions" | 9.193e8             | 7.602e8      | 3.025e9       | 1.277e10         | 8.587e8    | 
| "Number of X87"          | 1.048e4             | 1.718e4      | 4.168e4       | 8.762e4          | 7.262e4    | 
| "Number of SSE"          | 1.319e5             | 1.319e5      | 2.753e9       | 5.247e8          | 1.321e5    | 
| "Number of AVX"          | 7.088e8             | 6.882e8      | 0             | 5.736e9          | 8.197e8    | 
| "GFLOPS"                 | 21.12               | 4.607        | 1.155         | 11.04            | 3.267      | 
| "Efficiency"             | 71.36               | 15.56        | 3.904         | 37.31            | 11.04      | 
| "L1 Cache Loads"         | 3.413e8             | 1.51e6       | 2.235e6       | 5.097e9          | 3.434e7    | 
| "L1 Cache misses"        | 9.378e4             | 1.128e5      | 1.839e5       | 2.055e5          | 1.28e5     | 
| "Divider Active"         | 5.082e4             | 8.369e4      | 2.004e5       | 4.067e5          | 3.687e9    | 
| "Divider Active"         | 0.01025             | 0.0037       | 0.0022        | 0.0053           | 99.35      | 
| "No store buffer"        | 1.185e6             | 3.915e5      | 5.454e5       | 5.19e5           | 4.197e5    | 
| "No RS"                  | 2.581e8             | 1.996e9      | 8.018e9       | 4.06e9           | 3.41e9     | 
