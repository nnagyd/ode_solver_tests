using DifferentialEquations, DelimitedFiles, CPUTime, Plots

const numberOfParameters = 256
const numberOfRuns = 3
const numberOfTransient = 1024
const numberOfSaves = 32

#Parameter list
parameterList = collect(range(0.0,stop = 10.0,length = numberOfParameters))
q =  Vector{Float64}(undef,1)

#ODE
function valve!(dy,y,q,t)
    @inbounds begin
        dy[1] = y[2]
        dy[2] = -1.25*y[2]-(y[1]+10)+y[3]
        dy[3] = 20(q[1]-y[1]*sqrt(y[3]))
    end
    nothing
end

#event handling, when y[1] or y[2] changes sign
function condition(out,y,t,integrator)
    @inbounds begin
        out[1] = y[1] #valve is in contact with the seat of valve
        out[2] = y[2] #poincare section
    end
end

function local_max!(integrator, idx)
    @inbounds begin
        if idx == 1 #collision
            integrator.u[2] = -0.8*integrator.u[2]
        elseif idx == 2 #poincare section, end of iteration, start of new iteration
            terminate!(integrator)
        end
    end
end

function local_min!(integrator, idx)
    nothing
end

#defining ODE, Callback and ensemble simulation
cb = VectorContinuousCallback(
    condition,
    local_min!,2,
    affect_neg! = local_max!,
    rootfind = false, #does not work when true
    save_positions = (true,true),
    abstol=1e-6,reltol=1e-6)

#output data
global outputData = zeros(numberOfParameters,numberOfSaves*2);

#compile once
tSpan = [0.0,1e10]
y0 = [0.2,0.0,0.0]
q[1] = parameterList[1]
valveODE = ODEProblem(valve!,y0,tSpan,q)
res = solve(
    valveODE,
    DP5(),
    callback = cb,
    abstol = 1e-10,
    reltol = 1e-10,
    save_everystep = false,
    save_start = false,
    save_end = true,
    maxiters = 1e8,
    dense=false,
    dtmin = 1e-12)

#simulation
times =  Vector{Float64}(undef,numberOfRuns)
for runs in 1:numberOfRuns
    tStart = CPUtime_us()

    for par in 1:numberOfParameters
        #initial conditions
        tSpan = [0.0,1e10]
        y0 = [0.2,0.0,0.0]
        q[1] = parameterList[par]

        #transient simulation
        for i in 1:numberOfTransient
            valveODE = ODEProblem(valve!,y0,tSpan,q)
            res = solve(
                valveODE,
                DP5(),
                callback = cb,
                abstol = 1e-10,
                reltol = 1e-10,
                save_everystep = false,
                save_start = false,
                save_end = true,
                maxiters = 1e8,
                dense=false,
                dtmin = 1e-12)
            y0 = res.u[end]
            tSpan[1] = res.t[end]
        end

        #save values
        for i in 1:numberOfSaves
            valveODE = ODEProblem(valve!,y0,tSpan,q)
            res = solve(
                valveODE,
                DP5(),
                callback = cb,
                abstol = 1e-10,
                reltol = 1e-10,
                save_everystep = false,
                save_start = false,
                save_end = true,
                maxiters = 1e8,
                dense=false,
                dtmin = 1e-12)
            y0 = res.u[end]
            tSpan[1] = res.t[end]
            outputData[par,i] = res.u[end][1]
            outputData[par,i + 32] = res.u[end-2][1]
        end
    end
    tEnd = CPUtime_us()
    times[runs] = (tEnd-tStart)/(10^6)
    println("t = "*string(times[runs]))
 end

println("--------------------------------------------")
println("Time: "*string(times))
println("Parameter number: "*string(numberOfParameters))

plt = scatter(parameterList,outputData,marker = 1,color = :black,xlims =(0.,10.),ylims = (-0.2,10.),legend=false)

writedlm("data.csv", outputData, ',')
