using DifferentialEquations, DelimitedFiles, Plots, CPUTime, LoopVectorization

const numberOfRuns = 1
const numberOfParameters = 64
const unroll = 4
const numberOfTrajectories = Int32(numberOfParameters/unroll)

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
const P_A1_val = 1.5
const P_A2_val = 0
const f_1 = logRange(20.0,1_000.0,numberOfParameters)
const f_2 = 0

initialValues = Array{Float64,2}(undef,(numberOfTrajectories,2*unroll))
for i in 1:numberOfTrajectories
	for j in 0:unroll-1
		initialValues[i,1 + 2*j] = 1.0
	    initialValues[i,2 + 2*j] = 0.0
	end
end

#ODE settings
C = Vector{Float64}(undef,13*unroll) #C1-C13 -> ODE constants
y0 = initialValues[1,:] #inital conditions

#ODE system
function keller_miksis!(dy,y,C,τ)
	@avx for j = 1:unroll
		idf = j-1
		rx1 = 1/y[1 + 2idf]
        p = rx1 ^ C[10 + 13idf]

        s1,c1 = sincos(2*pi*τ)
        s2 = sin(2*pi*C[11 + 13idf]*τ+C[12 + 13idf])
		c2 = cos(2*pi*C[11 + 13idf]*τ+C[12 + 13idf])

        N1 = (C[13 + 13idf]+C[1 + 13idf]*y[2 + 2idf])*p-C[2 + 13idf]*(1+C[9 + 13idf]*y[ 2+ 2idf])-C[3 + 13idf]*rx1-C[4 + 13idf]*y[2 + 2idf]*rx1
		N2 = -1.5*(1.0-C[9+ 13idf]*y[2+ 2idf]*(1.0/3.0))*y[2 + 2idf]*y[2 + 2idf]
		N3 = -(C[5+ 13idf]*s1+C[6+ 13idf]*s2)*(1+C[9+ 13idf]*y[2 + 2idf])-y[1 + 2idf]*((C[7+ 13idf]*c1)+C[8+ 13idf]*c2)
		N = N1 + N2 + N3
        D = y[1+ 2idf]-C[9+ 13idf]*y[1+2idf]*y[2+2idf]+C[4+ 13idf]*C[9+ 13idf] #denominator
        rD = 1.0 / D

        #ODE system
        dy[1 + 2idf] = y[2 + 2idf]
        dy[2 + 2idf] = N * rD
	end
	nothing
end

#ensemble problem
function prob_func!(problem,i,repeat)
    @inbounds begin
		  for j = 0:unroll-1
			#calculating indexes
			problem.u0[1 + 2j] = initialValues[i,1 + 2j]
			problem.u0[2 + 2j] = initialValues[i,2 + 2j]

			#calculating physical parameters
			ω_1 = 2*pi*f_1[unroll*(i-1) + j + 1]*1000.0
			ω_2 = 0
			P_A1 = P_A1_val*1e5
			P_A2 = P_A2_val*1e5

			#calculating ODE constants
			tmp_1 = ((2*pi)/(R_E*ω_1))^2
			problem.p[1 + 13j] = (1-3*γ)/(ρ_L*c_L)*(P_inf - p_v + 2 * σ / R_E) * ((2*pi)/(R_E*ω_1))
			problem.p[2 + 13j] = (P_inf-p_v)/ρ_L*tmp_1
			problem.p[3 + 13j] = (2*σ)/(ρ_L*R_E)*tmp_1
			problem.p[4 + 13j] = 4*μ_L*2*pi/(ρ_L*R_E*R_E*ω_1)
			problem.p[5 + 13j] = P_A1/ρ_L*tmp_1
			problem.p[6 + 13j] = P_A2/ρ_L*tmp_1
			problem.p[7 + 13j] = R_E * ω_1 * P_A1/(ρ_L*c_L)*tmp_1
			problem.p[8 + 13j] = R_E * ω_1 * P_A2/(ρ_L*c_L)*tmp_1
			problem.p[9 + 13j] = R_E * ω_1 / (2*pi*c_L)
			problem.p[10 + 13j] = 3*γ
			problem.p[11 + 13j] = ω_2/ω_1
			problem.p[12 + 13j] = θ
			problem.p[13 + 13j] = (P_inf - p_v + 2 * σ / R_E)/ρ_L*tmp_1
		end
    end
    problem
end

# Compile once
tSpan = (0.0,1024.0)
prob = ODEProblem(keller_miksis!,y0,tSpan,C)
ensemble_prob = EnsembleProblem(prob,prob_func = prob_func!)
res = solve(
	ensemble_prob,
	DP5(),
	EnsembleSerial(),
	abstol = 1e-10,
	reltol = 1e-10,
	trajectories=numberOfTrajectories,
	save_everystep = false,
	save_start = false,
	save_end = true,
	dense = false,
	maxiters = 1e10,
	dtmin = 1e-10)
GC.gc()

#solving ODE 3x and measuring elapsed CPU time
times = Vector{Float64}(undef,numberOfRuns)
for runs = 1:numberOfRuns
    tStart = CPUtime_us()

    #transient
    global tSpan = (0.0,1024.0)
    global prob = ODEProblem(keller_miksis!,y0,tSpan,C)
    global ensemble_prob = EnsembleProblem(prob,prob_func = prob_func!)

    global res = solve(
        ensemble_prob,
        DP5(),
        EnsembleSerial(),
        abstol = 1e-10,
        reltol = 1e-10,
        trajectories= numberOfTrajectories,
        save_everystep = false,
        save_start = false,
        save_end = true,
        dense = false,
        maxiters = 1e10,
        dtmin = 1e-10)

    for i in 1:numberOfTrajectories
		for j in 1:2unroll
	        initialValues[i,j] = res[i].u[end][j]
		end
    end

    #save
    global tSpan = (0.0,64.0)
    global prob = ODEProblem(keller_miksis!,y0,tSpan,C)
    global ensemble_prob = EnsembleProblem(prob,prob_func = prob_func!)

    global res = solve(
        ensemble_prob,
        DP5(),
        EnsembleSerial(),
        abstol = 1e-10,
        reltol = 1e-10,
        trajectories= numberOfTrajectories,
        save_everystep = true,
        save_start = false,
        save_end = true,
        maxiters = 1e10,
        dense = false,
        dtmin = 1e-10
		)

    tEnd = CPUtime_us()
    times[runs] = (tEnd-tStart)/(10^6)
    println("t = "*string(times[runs]))
 end

println("----------------------")
println("Time: "*string(times))
println("Parameter number: "*string(numberOfParameters))

#finding global maxima
y_maxs = Vector{Float64}(undef,numberOfParameters)
for i = 1:numberOfTrajectories
	for j = 0:unroll-1
	    max = 0
	    for val in res[i].u
	        if max < val[1 + 2*j]
	            max = val[1 + 2*j]
	        end
	    end
	    y_maxs[1 + unroll*(i-1) + j] = max
	end
end

#plot for comparision
plotly()
scatter(f_1,y_maxs,xaxis = (:log,100.0:100.:1000),yaxis = 1:1:10,marker = (2,:black),legend = false,ylims = (1.,10.))
