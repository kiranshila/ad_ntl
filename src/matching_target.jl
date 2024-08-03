# Reading simulated data to build $\mathbf{T}$ from the paper

using StaticArrays, CSV

function read_ads_csv(csv_filename)
	file = CSV.File(csv_filename; comment = "#", types = [Float64, Float64, Float64, ComplexF64, ComplexF64, ComplexF64, ComplexF64, ComplexF64])
	freqs = file["freq"]
	n = length(freqs)
	# nfmin db to linear
	nfmin = @. 10^(file["nfmin"] / 10)
	# Scalar noise parameters
	rn = file["rn"]
	Sopt = file["sopt"]
	# Build S matrix of the target one freq at a time
	S = SMatrix{2, 2, ComplexF64}[]
	for i in 1:n
		push!(S, @SMatrix [
			file[i].s11 file[i].s12;
			file[i].s21 file[i].s22
		])
	end
	(freqs, S, nfmin, rn, Sopt)
end