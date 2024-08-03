using LinearAlgebra, StaticArrays, Statistics

include("$(@__DIR__)/../src/math.jl")

"""NTL Cascade"""
function cascade!(S_imn::AbstractArray{T}, zs, L, freqs, Z₀ = 50.0) where {T}
	δ = L / length(zs)
	tmp = similar(T)
	@inbounds @views for (i, freq) in enumerate(freqs)
		# Not really S Parmaters here, but we do everything inplace
		# So these are ABCD parameters until the end
		S_imn[i] .= T(I)
		ω = 2π * freq
		for z in zs
			abcd = tline2abcd(z, ω, δ)
			tmp .= S_imn[i]
			S_imn[i] .= tmp * abcd
		end
		S_imn[i] .= a2s(S_imn[i], Z₀)
	end
	S_imn
end

### Matching Problem Setup
#### Minimize the length of the NTL such that the magnitude^2 of S11 is less than some objective

function constraint!(cs, L::T1, zs::AbstractArray{T2}, freqs, Za, Zb) where {T1, T2}
	T = promote_type(T1, T2)
	# Big alloc here
	S_imn = [@MMatrix zeros(Complex{T}, 2, 2) for _ in 1:length(freqs)]
	# Solve the IMN
	cascade!(S_imn, zs, L, freqs, Za)
	# Compute the magnitude of the input reflection coefficient for the constraints
	Γl = (Zb - Za) / (Zb + Za)
	cs .= abs2.(Γin.(S_imn, Γl))
end

### Solving
using Optimization, OptimizationMOI, Ipopt, NLopt, CairoMakie

function solve_imn(N, freqs, Za, Zb, mΓ, Lhigh)
	Zlow = min(Za, Zb)    # Minimum impedance
	Zhigh = max(Za, Zb)  # Max impedance

	# Bounds and initial conditions
	u0 = [Lhigh, range(Za, Zb, N)...]
	lb = [0.0, fill(Zlow, N)...]
	ub = [Lhigh, fill(Zhigh, N)...]
	lcons = [fill(0.0, Nf)...]
	ucons = [fill(mΓ, Nf)...]

	f = OptimizationFunction((u, _, _...) -> u[1], AutoForwardDiff();
		cons = (res, x, _) -> constraint!(res, x[1], x[2:end], freqs, Za, Zb),
	)
	# IPopt
	opt = OptimizationMOI.MOI.OptimizerWithAttributes(Ipopt.Optimizer,
		"hessian_approximation" => "limited-memory")
	prob = OptimizationProblem(f, u0;
		lb = lb,
		ub = ub,
		lcons = lcons,
		ucons = ucons)
	solve(prob, opt)
end

Nf = 100
N = 20
Za = 50.0
Zb = 20.0
freqs = range(10, 100.0, Nf) .* 1e9
sol = solve_imn(N, freqs, Za, Zb, 0.0001, 30e-3)

gamma_ins = zeros(Nf)
constraint!(gamma_ins, sol.u[1], sol.u[2:end], freqs, Za, Zb)
lines(freqs, @. 10 * log10(gamma_ins))
plot(range(0, sol.u[1], N), sol.u[2:end])
