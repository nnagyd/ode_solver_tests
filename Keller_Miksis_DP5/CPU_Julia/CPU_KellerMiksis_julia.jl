using DifferentialEquations, DelimitedFiles, Plots, CPUTime

const numberOfRuns = 3
const numberOfParameters = 16 #number of different parameters for each control variable

function logRange(startVal,endVal,intervals)
    intervalDelta = endVal/startVal
    λ = intervalDelta^(1/(intervals-1))
    newVal = startVal
    vals = [newVal]
    for i in 1:intervals-1
        newVal *= λ
        push!(vals,newVal)
    end
    vals
end

#physical parameters
const ρ_L = 9.970639504998557e+02
const P_inf = 1.0e+5
const p_v = 3.166775638952003e+03
const σ = 0.071977583160056
const R_E = 10/1e6
const γ = 1.4
const c_L = 1.497251785455527e+03
const μ_L = 8.902125058209557e-04
const θ = 0

#outer parameters are P_A1 and P_A2
P_A1_val = 1.5 #1.1
P_A2_val = 0 #1.2
f_1 = logRange(20.0,1_000.0,numberOfParameters)
f_2 = 0

initialValues = Array{Float64,2}(undef,(numberOfParameters,2))
for i in 1:numberOfParameters
    initialValues[i,1] = 1.0
    initialValues[i,2] = 0.0
end

#ODE settings
C = Vector{Float64}(undef,13) #C1-C13 -> ODE constants
y0 = [1.0,0.0] #inital conditions0

#ODE system
function keller_miksis!(dy,y,C,τ)
    @inbounds begin
        #parts of nominator
        p1 = (C[13]+C[1]*y[2])*(1/y[1])^C[10]-C[2]*(1+C[9]*y[2])-C[3]/y[1]-C[4]*y[2]/y[1]
        p2 = -(1-C[9]*y[2]/3.0)*3/2*y[2]*y[2]-(C[5]*sin(2*pi*τ)+C[6]*sin(2*pi*C[11]*τ+C[12]))*(1+C[9]*y[2])
        p3 = -y[1]*((C[7]*cos(2*pi*τ))+C[8]*cos(2*pi*C[11]*τ+C[12]))

        N_KM = p1+p2+p3 #nominator
        D_KM = y[1]-C[9]*y[1]*y[2]+C[4]*C[9] #denominator

        #ODE system
        dy[1] = y[2]
        dy[2] = N_KM/D_KM
    end
    nothing
end

#ensemble problem
function prob_func!(problem,i,repeat)
    @inbounds begin
        println(i)
        #calculating indexes
        problem.u0[1] = initialValues[i,1]
        problem.u0[2] = initialValues[i,2]

        #calculating physical parameters
        ω_1 = 2*pi*f_1[i]*1000.0
        ω_2 = 0
        P_A1 = P_A1_val*1e5
        P_A2 = P_A2_val*1e5

        #calculating ODE constants
        tmp_1 = ((2*pi)/(R_E*ω_1))^2
        problem.p[1] = (1-3*γ)/(ρ_L*c_L)*(P_inf - p_v + 2 * σ / R_E) * ((2*pi)/(R_E*ω_1))
        problem.p[2] = (P_inf-p_v)/ρ_L*tmp_1
        problem.p[3] = (2*σ)/(ρ_L*R_E)*tmp_1
        problem.p[4] = 4*μ_L*2*pi/(ρ_L*R_E*R_E*ω_1)
        problem.p[5] = P_A1/ρ_L*tmp_1
        problem.p[6] = P_A2/ρ_L*tmp_1
        problem.p[7] = R_E * ω_1 * P_A1/(ρ_L*c_L)*tmp_1
        problem.p[8] = R_E * ω_1 * P_A2/(ρ_L*c_L)*tmp_1
        problem.p[9] = R_E * ω_1 / (2*pi*c_L)
        problem.p[10] = 3*γ
        problem.p[11] = ω_2/ω_1
        problem.p[12] = θ
        problem.p[13] = (P_inf - p_v + 2 * σ / R_E)/ρ_L*tmp_1
    end
    problem
end

#solving ODE 3x and measuring elapsed CPU time
times = Vector{Float64}(undef,numberOfRuns)
for runs in 1:numberOfRuns
    tStart = CPUtime_us()

    #transient
    tSpan = (0.0,1024.0)
    prob = ODEProblem(keller_miksis!,y0,tSpan,C)
    ensemble_prob = EnsembleProblem(prob,prob_func = prob_func!)

    global res = solve(
        ensemble_prob,
        DP5(),
        EnsembleSerial(),
        abstol = 1e-10,
        reltol = 1e-10,
        trajectories= numberOfParameters,
        save_everystep = false,
        save_start = false,
        save_end = true,
        dense = false,
        maxiters = 1e10,
        dtmin = 1e-10)

    for i in 1:numberOfParameters
        initialValues[i,1] = res[i].u[end][1]
        initialValues[i,2] = res[i].u[end][2]
    end

    #save
    tSpan = (0.0,64.0)
    prob = ODEProblem(keller_miksis!,y0,tSpan,C)
    ensemble_prob = EnsembleProblem(prob,prob_func = prob_func!)

    global res = solve(
        ensemble_prob,
        DP5(),
        EnsembleSerial(),
        abstol = 1e-10,
        reltol = 1e-10,
        trajectories= numberOfParameters,
        save_everystep = true,
        save_start = false,
        save_end = true,
        maxiters = 1e10,
        dense = false,
        dtmin = 1e-10)

    tEnd = CPUtime_us()
    times[runs] = (tEnd-tStart)/(10^6)
    println("t = "*string(times[runs]))
 end

println("----------------------")
println("Time: "*string(times))
println("Parameter number: "*string(numberOfParameters))

outputData = zeros(numberOfParameters,4)
for i in 1:numberOfParameters
    outputData[i,1] = i
    outputData[i,2] = f_1[i]
    outputData[i,3] = res[i].u[end][1]
    outputData[i,4] = res[i].u[end][2]
end

writedlm("keller_miksis_endvalues.csv", outputData, ',')

y_maxs = Vector{Float64}(undef,numberOfParameters)
for i in 1:numberOfParameters
    max = 0
    for val in res[i].u
        if max < val[1]
            max = val[1]
        end
    end
    y_maxs[i] = max
end

scatter(f_1,y_maxs,xaxis = (:log,100.0:100.:1000),yaxis = 1:1:10,marker = (2,:black),legend = false,ylims = (1.,10.))
