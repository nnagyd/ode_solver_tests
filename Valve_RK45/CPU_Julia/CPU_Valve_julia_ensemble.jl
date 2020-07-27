using DifferentialEquations, DelimitedFiles, Plots, CPUTime

const numberOfParameters = 256
const numberOfRuns = 3

function valve!(dy,y,q,t)
    @inbounds begin
        dy[1] = y[2]
        dy[2] = -1.25*y[2]-(y[1]+10)+y[3]
        dy[3] = 20(q[1]-y[1]*sqrt(y[3]))
    end
    nothing
end

tSpan = (0.0,1e10)
parameterList = collect(range(0.2,stop = 10.0,length = numberOfParameters))
q =  Vector{Float64}(undef,4) #q[1] --> q parameter q[2] --> number of iterations q[3] --> number of equation, q[4] last extremum
y0 = [0.2,0.0,0.0]
global outputData = ones(numberOfParameters,65)*-1 #initialize output array with -1s

#ensemble problem
function parameter_change!(problem,i,repeat)
    println(i)
    @inbounds begin
        problem.p[1] = parameterList[i] #parameter
        problem.p[3] = i #number
        global outputData[i,1] = parameterList[i] #saving to global matrix
    end
    problem
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
            integrator.p[4] = integrator.u[1] #save this local minima for convergence check
            iterationNumber = integrator.p[2]
            if iterationNumber > 1024
                i = Int64(integrator.p[3])
                j = Int64(iterationNumber - 1024);
                outputData[i,Int64(j+33)] = integrator.u[1]
            end
        elseif idx == 2 #poincare section, end of iteration, start of new iteration
            integrator.p[2] = integrator.p[2]+1
            iterationNumber = integrator.p[2]
            i = Int64(integrator.p[3])

            if iterationNumber > 1024
                j = Int64(iterationNumber - 1024);
                outputData[i,j+1] = integrator.u[1]
            end

            #max iterations reached
            if iterationNumber >= 1024+32
                integrator.t = 1e10
            end

            # convergence detected at u[1] != 0
            if abs(integrator.p[4] - integrator.u[1]) < 1e-9 && abs(integrator.u[1])>1e-9
                outputData[i,2] = integrator.u[1]
                integrator.t = 1e10
            end

            integrator.p[4] = integrator.u[1]
        end
    end
end

function local_min!(integrator, idx)
    @inbounds begin
        if idx == 2 #detection of local minima, save it when 1024 iterations are reached
            iterationNumber = integrator.p[2]
            integrator.p[4] = integrator.u[1] #save this local minima for convergence check
            if iterationNumber > 1024
                i = Int64(integrator.p[3])
                j = Int64(iterationNumber - 1024);
                outputData[i,Int64(j+33)] = integrator.u[1]
            end
        end
    end
end

#defining ODE, Callback and ensemble simulation
valveODE = ODEProblem(valve!,y0,tSpan,q)
cb = VectorContinuousCallback(condition,local_min!,2,affect_neg! = local_max!)
ensemble_prob = EnsembleProblem(valveODE,prob_func = parameter_change!)


times =  Vector{Float64}(undef,numberOfRuns)
for runs in 1:numberOfRuns
    tStart = CPUtime_us()
    global res = solve(
        ensemble_prob,
        DP5(),
        EnsembleSerial(),
        callback = cb,
        abstol = 1e-10,
        reltol = 1e-10,
        trajectories= numberOfParameters,
        save_everystep = false,
        save_start = false,
        save_end = true,
        maxiters = 1e8,
        dtmin = 1e-12)
    tEnd = CPUtime_us()
    times[runs] = (tEnd-tStart)/(10^6)
    println("t = "*string(times[runs]))
 end

println("----------------------")
println("Time: "*string(times))
println("Parameter number: "*string(numberOfParameters))

plt = scatter(marker = 2,legend = false,xlims =(0.,10.),ylims = (-0.2,10.))
for j in 2:65
    scatter!(plt,outputData[:,1],outputData[:,j],marker = 1,color = :black)
end

writedlm("data.csv", outputData, ',')
