using DifferentialEquations, CPUTime, Statistics

#settings
const unroll = 128
const numberOfParameters = 46080
const numberOfRuns = 3

#parameters and initial conditions
parameterList = range(0.0, stop = 21.0, length = numberOfParameters)
p = zeros(unroll,1)
tspan = (0.0, 10.0)
u0 = ones(3*unroll,1)*10

#ODE
function lorenz!(du, u, p, t)
  for j in 1:128
    index = 3(j-1)+1
    du[index] = 10.0 * (u[index+1] - u[index])
    du[index+1] = p[j] * u[index] - u[index+1] - u[index] * u[index+2]
    du[index+2] = u[index] * u[index+1] - 2.66666666 * u[index+2]
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
    save_end = true,
    adaptive = false,
    dt = 0.01,
    dense = false
    )

    tEnd = CPUtime_us()
    times[runs] = (tEnd-tStart)/(10^6)
end

println("End of test")
println("Parameter number: "*string(numberOfParameters))
println("Unroll: "*string(unroll))
println("Time: "*string(times[2:numberOfRuns]))
println("Avg: "*string(mean(times[2:numberOfRuns])))
