using MKL # Acceleration on Intel machines
using LinearAlgebra, StaticArrays, Statistics

include("$(@__DIR__)/../src/math.jl")
include("$(@__DIR__)/../src/interp_rlgc.jl")
include("$(@__DIR__)/../src/matching_target.jl")
include("$(@__DIR__)/../src/dxf.jl")

# Same code as the microstrip example, this time with the WBLNA SS cross-section
# Don't convert the output to S, as we need to cascade with the matching target, though
function cascade!(ABCD_imn::AbstractArray{T}, widths, L, freqs) where {T}
	δ = L / length(widths)
	tmp = similar(T)
	@inbounds @views for (i, freq) in enumerate(freqs)
		ABCD_imn[i] .= T(I)
		ω = 2π * freq
		for w in widths
			rlgc = WBLNA_SS_MODEL(freq, w)
			abcd = rlgc2abcd(rlgc, ω, δ)
			tmp .= ABCD_imn[i]
			ABCD_imn[i] .= tmp * abcd
		end
	end
	ABCD_imn
end

# Compute the behavior of the IMN connector to the matching target
# This manifests as noise and input return loss
function evaluate(L::T1, widths::AbstractVector{T2}, freqs, ABCD_T, NFmin_T, Rn_T, Sopt_T) where {T1, T2}
	T = promote_type(T1, T2)

	# Compute the IMN
	ABCD_imn = [@MMatrix zeros(Complex{T}, 2, 2) for _ in 1:length(freqs)]
	cascade!(ABCD_imn, widths, L, freqs)

	# Compute the S11 of the cascade with the target amplifier
	ABCD_casc = @. ABCD_imn * ABCD_T
	S_casc = @. SMatrix(a2s(ABCD_casc))
	S11 = [abs2(S[1, 1]) for S in S_casc]

	# Compute the noise due to mismatch
	S_imn = @. SMatrix(a2s(ABCD_imn))
	# Assume a perfect match at the IMN input, compute the impedance presented to the target's input
	gamma_out = Γout.(S_imn, 0.0)
	# Compute the input noise of the amplifier via the noise parameters and the source impedance
	NF_amp = noise_figure.(NFmin_T, Rn_T, Sopt_T, gamma_out)
	# Convert the noise figure to noise temperature
	Te_amp = noise_temperature.(NF_amp)

	# Compute the equivalent noise of the (passive) IMN, assuming a physical temperature of 290
	# and convert to noise temperature
	NF_imn = noise_figure.(S_imn, 0.0, 290)
	Te_imn = noise_temperature.(NF_imn)
	# Convert NF_min to Te_min
	Te_amp_min = noise_temperature.(NFmin_T)

	# Add the two noises to get the input noise and the offset tmin
	S11, Te_amp + Te_imn, Te_amp_min + Te_imn
end

### Optimization Setup
#### Minimize the mean noise of the amplfiier such that the magnitude^2 of S11 is less than some objective

function objective(L, widths, freqs, ABCD_T, NFmin_T, Rn_T, Sopt_T)
	# Compute the noise
	_, t, _ = evaluate(L, widths, freqs, ABCD_T, NFmin_T, Rn_T, Sopt_T)
	# Find the band average
	mean(t)
end

function constraint!(cs, L, widths, freqs, ABCD_T, NFmin_T, Rn_T, Sopt_T)
	# Compute the input return loss (already abs^2-ed)
	S11, _, _ = evaluate(L, widths, freqs, ABCD_T, NFmin_T, Rn_T, Sopt_T)
	# Create a constraint for every point in frequency in dB
	cs .= @. 10 * log10(S11)
end

### Solving
using Optimization, OptimizationMOI, Ipopt, CairoMakie
using FiniteDiff # To compare

function solve_imn(N, S11_dB::Float64; kwargs...)
	# Extract data for T
	freqs, S_T, NFmin_T, Rn_T, Sopt_T = read_ads_csv("$(@__DIR__)/../matching_target_data/LNA_Noise_S.csv")
	# Precompute ABCD for the target, as we are going to cascade it with ABCD on the IMN
	ABCD_T = @. SMatrix(s2a(S_T))
	Nf = length(freqs)

	# Geometric constraints
	Wlow = 0.1e-3
	Whigh = 15e-3

	# Bounds and initial conditions
	u0 = [range(Whigh, Wlow, N)..., 100e-3]
	#u0 = [((Whigh-Wlow).*rand(N) .+ Wlow)..., 100e-3] # Random starting position
	lb = [fill(Wlow, N)..., 50e-3]
	ub = [fill(Whigh, N)..., 120e-3]
	# Constraint of -Inf < |S11| <= S11_dB
	lcons = [fill(-Inf, Nf)...]
	ucons = [fill(S11_dB, Nf)...]

	# Create version of the objective and constraints curried with locals
	curried_objective = (u, _, _...) -> objective(u[N+1], u[1:N], freqs, ABCD_T, NFmin_T, Rn_T, Sopt_T)
	curried_constraint = (cs, u, _) -> constraint!(cs, u[N+1], u[1:N], freqs, ABCD_T, NFmin_T, Rn_T, Sopt_T)

	# Setup the optimization problem
	# f = OptimizationFunction(curried_objective, AutoFiniteDiff(); cons = curried_constraint)
	f = OptimizationFunction(curried_objective, AutoForwardDiff(); cons = curried_constraint)
	opt = OptimizationMOI.MOI.OptimizerWithAttributes(Ipopt.Optimizer, "hessian_approximation" => "limited-memory")
	prob = OptimizationProblem(f, u0; lb = lb, ub = ub, lcons = lcons, ucons = ucons)

	# And solve
	sol = solve(prob, opt; kwargs...)
	L = sol.u[end]
	widths = sol.u[1:N]

	# Return the things we care about (the design and behavior)
	S11, Te, Tmin = evaluate(L, widths, freqs, ABCD_T, NFmin_T, Rn_T, Sopt_T)
	(widths, L, S11, Te, Tmin, freqs)
end

N = 100
Γmax = -9.3
widths, L, S11, Te, Tmin, freqs = solve_imn(N, Γmax; maxtime = 60)

### Visualization and Exporting Design

begin
	f = Figure(size = (800, 600))

	# Plot the shape
	ax = Axis(f[2, 1:2], aspect = DataAspect(), title = "NTL Profile", limits = (0, 100, -15, 15))
	stairs!(ax, range(0, L, N) .* 1000, -widths ./ 2 .* 1000, color = :black, step = :center)
	stairs!(ax, range(0, L, N) .* 1000, widths ./ 2 .* 1000, color = :black, step = :center)

	# Plot the performance
	ax = Axis(f[1, 1], ylabel = "S11 (dB)", xlabel = "Freq (GHz)", title = "Return Loss", limits = (0.7, 2, -20, 0))
	lines!(ax, freqs ./ 1e9, @. 10 * log10(S11))
	μT = mean(Te)
	ax = Axis(f[1, 2], ylabel = "Noise (K)", xlabel = "Freq (GHz)", title = "Noise: $(round(μT,digits=2)) K Average", limits = (0.7, 2, 0, 20))
	lines!(ax, freqs ./ 1e9, Te, label = "Tamp")
	lines!(ax, freqs ./ 1e9, Tmin, label = "Tmin")
	axislegend(ax, position = :rt)

	f
end

# Generate geometry file
create_dxf(L, widths, "$(@__DIR__)/../shape_results/N$(N)_$(Γmax).dxf")