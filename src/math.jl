"""Speed of light in a vacuum"""
const c₀ = 299_792_458

"""Vacuum permeability"""
const μ₀ = 1.25663706127e-6

"""Vacuum permitivity"""
const ε₀ = 8.854e-12

"""Imepdance of free space"""
const η₀ = sqrt(μ₀ / ε₀)

"""ABCD to S Parameters"""
function a2s(ABCD, Z₀ = 50.0)
	denom = ABCD[1] + ABCD[3] / Z₀ + ABCD[2] * Z₀ + ABCD[4]
	S11 = (ABCD[1] + ABCD[3] / Z₀ - ABCD[2] * Z₀ - ABCD[4]) / denom
	S21 = 2 * (ABCD[1] * ABCD[4] - ABCD[3] * ABCD[2]) / denom
	S12 = 2 / denom
	S22 = (-ABCD[1] + ABCD[3] / Z₀ - ABCD[2] * Z₀ + ABCD[4]) / denom
	@SMatrix [S11 S21; S12 S22]
end

"""S to ABCD Parameters"""
function s2a(S::AbstractMatrix, Z₀ = 50.0)
	A = ((1 + S[1, 1]) * (1 - S[2, 2]) + S[1, 2] * S[2, 1]) / (2 * S[2, 1])
	B = Z₀ * ((1 + S[1, 1]) * (1 + S[2, 2]) - S[1, 2] * S[2, 1]) / (2 * S[2, 1])
	C = 1 / Z₀ * ((1 - S[1, 1]) * (1 - S[2, 2]) - S[1, 2] * S[2, 1]) / (2 * S[2, 1])
	D = ((1 - S[1, 1]) * (1 + S[2, 2]) + S[1, 2] * S[2, 1]) / (2 * S[2, 1])
	@SMatrix [A B; C D]
end

"""Input reflection coefficient"""
Γin(S, Γl) = S[1, 1] + (S[2, 1] * S[1, 2] * Γl) / (1 - S[2, 2] * Γl)

"""Output reflection coefficient"""
Γout(S, Γs) = S[2, 2] + (S[1, 2] * S[2, 1] * Γs) / (1 - S[1, 1] * Γs)

"""ABCD matrix from transmission line parameters"""
function tline2abcd(Z₀, ω, l; α = 0.0, νₚ = c₀)
	β = ω / νₚ
	γ = α + β * im
	sh = sinh(γ * l)
	ch = cosh(γ * l)
	@SMatrix [ch sh*Z₀; sh/Z₀ ch]
end

struct RLGC{T <: Real}
	R::T
	L::T
	G::T
	C::T
end

characteristic_impedance(rlgc::RLGC, ω) = sqrt((rlgc.R + im * ω * rlgc.L) / (rlgc.G + im * ω * rlgc.C))
propagation_constant(rlgc::RLGC, ω) = sqrt((rlgc.R + im * ω * rlgc.L) * (rlgc.G + im * ω * rlgc.C))

"""ABCD matrix from RLGC parameters"""
function rlgc2abcd(rlgc::RLGC, ω, l)
	γ = propagation_constant(rlgc, ω)
	z₀ = characteristic_impedance(rlgc, ω)
	sh = sinh(γ * l)
	ch = cosh(γ * l)
	@SMatrix [ch sh*z₀; sh/z₀ ch]
end

function available_gain(S, Γs)
	# Pozar eqn. 5.88
	pavn = abs2(S[2, 1]) * (1 - abs2(Γs))
	pavs = abs2(1 - S[1, 1] * Γs) * (1 - abs2(Γout(S, Γs)))
	pavn / pavs
end

"""IEEE Standard Reference Temperature in Kelvin"""
const T₀ = 290

"""Noise figure from the four noise paramters and a source reflection Γs"""
function noise_figure(nfmin, rn, Γopt, Γs; z0 = 50.0)
	nfmin + (4 * rn * (abs2(Γs - Γopt))) / (z0 * (1 - abs2(Γs)) * abs2(1 + Γopt))
end

"""Noise figure of a passive two port network (linear)"""
function noise_figure(S, Γs, T)
	ga = available_gain(S, Γs)
	1 + ((1 - ga) / ga) * (T / T₀)
end

"""Convert noise figure (linear) to noise temperature"""
noise_temperature(nf) = T₀ * (nf - 1)
