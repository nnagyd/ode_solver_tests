using DifferentialEquations, CPUTime, Statistics

#settings
const rollOut = 128
const numberOfParameters = 46080
const numberOfRuns = 3

#parameters and initial conditions
parameterList = collect(range(0.0, stop = 21.0, length = numberOfParameters))
p = Vector{Float64}(undef,rollOut)
tspan = (0.0, 10.0)
u0 = Vector{Float64}(undef,3rollOut)
for i in 1:3rollOut
  u0[i] = 10.
end

#ODE
function lorenz!(du, u, p, t)
  @inbounds begin
    @simd for j in 1:rollOut
      index = 3(j-1)+1
      du[index] = 10.0 * (u[index+1] - u[index])
      du[index+1] = p[j] * u[index] - u[index+1] - u[index] * u[index+2]
      du[index+2] = u[index] * u[index+1] - 2.66666666 * u[index+2]
    end
  end
end


#simulation
times = Vector{Float64}(undef,numberOfRuns)
for runs in 1:numberOfRuns
  tStart = CPUtime_us()

  for j in 1:rollOut:numberOfParameters
    #defining ODE
    p = parameterList[j:j+rollOut-1]
    prob = ODEProblem(lorenz!, u0, tspan, p)

    solve(
      prob,
      RK4(),
      save_everystep = false,
      save_start = false,
      save_end = true,
      adaptive = false,
      dt = 0.01,
      )
  end

  tEnd = CPUtime_us()
  times[runs] = (tEnd-tStart)/(10^6)
end

println("End of test")
println("Parameter number: "*string(numberOfParameters))
println("RollOut: "*string(rollOut))
println("Time: "*string(times[2:numberOfRuns]))
println("Avg: "*string(mean(times[2:numberOfRuns])))
