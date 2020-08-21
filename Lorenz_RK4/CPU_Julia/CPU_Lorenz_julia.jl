using DifferentialEquations, CPUTime, Statistics, MuladdMacro

#settings
const unroll = 128
const numberOfParameters = 46080
const numberOfRuns = 3
nocheck(dt,u,p,t) = false

#parameters and initial conditions
parameterList = collect(range(0.0, stop = 21.0, length = numberOfParameters))
p = Vector{Float64}(undef,unroll)
tspan = (0.0, 10.0)
u0 = Vector{Float64}(undef,3unroll)
for i in 1:3unroll
  u0[i] = 10.
end

#ODE
function lorenz!(du, u, p, t)
    @muladd @inbounds for j in 1:unroll
      du[j] = 10.0 * (u[j+unroll] - u[j])
      du[j+unroll] = p[j] * u[j] - u[j+unroll] - u[j] * u[j+2*unroll]
      du[j+2*unroll] = u[j] * u[j+unroll] - 2.66666666 * u[j+2*unroll]
   end
end

#compile
p = parameterList[1:1+unroll-1]
prob = ODEProblem(lorenz!, u0, tspan, p)

solve(
  prob,
  RK4(),
  save_everystep = false,
  save_start = false,
  save_end = true,
  adaptive = false,
  dt = 0.01,
  unstable_check = nocheck
  )
GC.gc()

#simulation
times = Vector{Float64}(undef,numberOfRuns)
for runs in 1:numberOfRuns
  tStart = CPUtime_us()

  for j in 1:unroll:numberOfParameters
    #defining ODE
    p = parameterList[j:j+unroll-1]
    prob = ODEProblem(lorenz!, u0, tspan, p)

    solve(
      prob,
      RK4(),
      save_everystep = false,
      save_start = false,
      save_end = true,
      adaptive = false,
      dt = 0.01,
      unstable_check = nocheck
      )
  end

  tEnd = CPUtime_us()
  times[runs] = (tEnd-tStart)/(10^6)
end

println("End of test")
println("Parameter number: "*string(numberOfParameters))
println("Unroll: "*string(unroll))
println("Time: "*string(times[2:numberOfRuns]))
println("Avg: "*string(mean(times[2:numberOfRuns])))
