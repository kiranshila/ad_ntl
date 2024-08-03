using LinearAlgebra, StaticArrays, Statistics

include("$(@__DIR__)/../src/interp_rlgc.jl")
include("$(@__DIR__)/../src/math.jl")

function cascade!(S_imn::AbstractArray{T}, widths, L, freqs, Z₀ = 50.0) where {T}
	δ = L / length(widths)
	tmp = similar(T)
	@inbounds @views for (i, freq) in enumerate(freqs)
		S_imn[i] .= T(I)
		ω = 2π * freq
		for w in widths
			rlgc = RO5880_MODEL(freq, w)
			abcd = rlgc2abcd(rlgc, ω, δ)
			tmp .= S_imn[i]
			S_imn[i] .= tmp * abcd
		end
		S_imn[i] .= a2s(S_imn[i], Z₀)
	end
	S_imn
end

### Matching Problem Setup
#### Minimize the length of the NTL such that the magnitude^2 of S11 is less than some objective

objective(L) = L

function constraint!(cs, widths::AbstractArray{T1}, L::T2, freqs, Za, Zb) where {T1, T2}
	T = promote_type(T1, T2)
	S_imn = [@MMatrix zeros(Complex{T}, 2, 2) for _ in 1:length(freqs)]
	cascade!(S_imn, widths, L, freqs, Za)
	Γl = (Zb - Za) / (Zb + Za)
	cs .= abs2.(Γin.(S_imn, Γl))
end

### Solving
using Optimization, OptimizationMOI, Ipopt, NLopt, CairoMakie

function solve_imn(N, freqs, Za, Zb, mΓ)
	# We want to keep W/H greater than 1
	Wlow = 0.252e-3
	Whigh = 3e-3

	# Set the initial length to lambda /2
	λ_max = c₀ / sqrt(2.2) / minimum(freqs)
	L = λ_max / 2

	# Bounds and initial conditions
	u0 = [fill(1e-3, N)..., L * 2]
	lb = [fill(Wlow, N)..., L]
	ub = [fill(Whigh, N)..., Inf]
	lcons = [fill(0.0, Nf)...]
	ucons = [fill(mΓ, Nf)...]

	f = OptimizationFunction((u, _, _...) -> objective(u[N+1]), AutoForwardDiff();
		cons = (res, x, _) -> constraint!(res, x[1:N], x[N+1], freqs, Za, Zb),
	)
	# IPopt
	opt = OptimizationMOI.MOI.OptimizerWithAttributes(Ipopt.Optimizer, "hessian_approximation" => "limited-memory")
	prob = OptimizationProblem(f, u0;
		lb = lb,
		ub = ub,
		lcons = lcons,
		ucons = ucons)
	solve(prob, opt)
end

Nf = 71
N = 20
Za = 50.0
Zb = 20.0
freqs = range(1, 10.0, Nf) .* 1e9
sol = solve_imn(N, freqs, Za, Zb, 0.001)

gamma_ins = zeros(Nf)
constraint!(gamma_ins, sol.u[1:N], sol.u[N+1], freqs, Za, Zb)
lines(freqs, @. 10 * log10(gamma_ins))
plot(range(0, sol.u[N+1], N), sol.u[1:N])
