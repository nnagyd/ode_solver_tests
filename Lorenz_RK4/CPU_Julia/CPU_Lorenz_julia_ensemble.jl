using DifferentialEquations, CPUTime, Statistics, LoopVectorization

#settings
const unroll = 128
const numberOfParameters = 46080
const numberOfRuns = 3
nocheck(dt,u,p,t) = false

#parameters and initial conditions
parameterList = range(0.0, stop = 21.0, length = numberOfParameters)
p = zeros(unroll,1)
tspan = (0.0, 10.0)
u0 = ones(3*unroll,1)*10

#ODE
function lorenz!(du, u, p, t)
    @avx for j in 1:unroll
      du[j] = 10.0 * (u[j+unroll] - u[j])
      du[j+unroll] = p[j] * u[j] - u[j+unroll] - u[j] * u[j+2*unroll]
      du[j+2*unroll] = u[j] * u[j+unroll] - 2.66666666 * u[j+2*unroll]
   end
end

#parameter change function
function parameterChange!(prob, i, repeat)
  index = Int64(unroll*i)
  for j in 1:unroll
    prob.p[j] = parameterList[index-j+1]
  end
  prob
end

#defining ode problem
lorenzProblem = ODEProblem(lorenz!, u0, tspan, p)
ensemble_prob = EnsembleProblem(lorenzProblem, prob_func = parameterChange!)

#compile
solve(
  ensemble_prob,
  RK4(),
  EnsembleSerial(),
  trajectories = numberOfParameters/unroll,
  save_everystep = false,
  save_start = false,
  adaptive = false,
  dt = 0.01,
  dense = false,
  unstable_check = nocheck,
  )
GC.gc()

#simulation
times = Vector{Float64}(undef,numberOfRuns)
for runs in 1:numberOfRuns
  tStart = CPUtime_us()
  solve(
    ensemble_prob,
    RK4(),
    EnsembleSerial(),
    trajectories = numberOfParameters/unroll,
    save_everystep = false,
    save_start = false,
    adaptive = false,
    dt = 0.01,
    dense = false,
    unstable_check = nocheck,
    )

    tEnd = CPUtime_us()
    times[runs] = (tEnd-tStart)/(10^6)
end

println("End of test")
println("Parameter number: "*string(numberOfParameters))
println("Unroll: "*string(unroll))
println("Time: "*string(times[2:numberOfRuns]))
println("Avg: "*string(mean(times[2:numberOfRuns])))
