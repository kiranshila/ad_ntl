using Interpolations, CSV

# We're going to interpolate the same way every time
interp(ys, xs) = extrapolate(interpolate(xs, ys, Gridded(Linear())), Line())

"""Generate a causal transmission line model based on simulated RLGC data"""
function build_interp_rlgc(filename)
	# Read the data
	rlgc_data = CSV.File(filename; comment = "#")
	# Determine the data ranges and scale
	ws = unique(rlgc_data["w_mm"]) .* 1e-3
	fs = unique(rlgc_data["freq_ghz"]) .* 1e9
	data_points = (fs, ws)
	data_shape = (length(fs), length(ws))
	# Extract and scale
	rlgc_r = reshape(rlgc_data["r_ohm"], data_shape)
	rlgc_l = if :l_uH in rlgc_data.names
		reshape(rlgc_data["l_uH"], data_shape) .* 1e-6
	elseif :l_nH in rlgc_data.names
		reshape(rlgc_data["l_nH"], data_shape) .* 1e-9
	else
		throw("Bad L scaling factor")
	end
	rlgc_g = reshape(rlgc_data["g_uS"], data_shape) .* 1e-6
	rlgc_c = reshape(rlgc_data["c_pF"], data_shape) .* 1e-12
	# Create the interpolations
	r_model = interp(rlgc_r, data_points)
	l_model = interp(rlgc_l, data_points)
	g_model = interp(rlgc_g, data_points)
	c_model = interp(rlgc_c, data_points)
	# And return the function that builds the RLGC parameters
	(f, w) -> begin
		r = r_model(f, w)
		l = l_model(f, w)
		g = g_model(f, w)
		c = c_model(f, w)
		RLGC(r, l, g, c)
	end
end

const RO5880_MS_MODEL = build_interp_rlgc("$(@__DIR__)/../rlgc_data/ro5880_10mil_rlgc.csv")
const WBLNA_SS_MODEL = build_interp_rlgc("$(@__DIR__)/../rlgc_data/wblna_7x10.csv")
