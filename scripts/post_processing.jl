using Measurements, CSV, DataFrames, SavitzkyGolay, CairoMakie, Statistics

# Mismatch and available gain correction factor
correction(Γₛ, Γₐ) = (1 - abs2(Γₛ)) / abs2(1 - Γₛ * Γₐ)

# Curve Fit TCold Model, interpreting errors as 3σ
tcold(f_ghz) = 78.0824 + 1.2425 * f_ghz - 0.1166 * f_ghz * f_ghz
tcold_err(f_ghz) = (0.0776 + 0.0568 * f_ghz - 0.0048 * f_ghz * f_ghz) / 3
tcold_meas(f_ghz) = tcold(f_ghz) ± tcold_err(f_ghz)

# Utility
dbm2w(dbm) = 10^(dbm / 10) / 1e3

# Compute the mean and and uncertainty of the SA power measurements
function col_meas(raw, col; N=5)
    cols = ["$(col)_$(i)" for i in 0:(N-1)]
    block = Matrix(select(raw, cols))
    means = mean(block, dims=2)
    vec(@. dbm2w(means ± 0.02 / sqrt(N)))
end

# Curve Fit PNA Uncertainty Model, interpreting errors as 3σ
phase_err(mag) = (0.25 / mag + 0.50) / 3
mag_err(mag) = (0.004 + 0.006 * (mag + mag^2)) / 3
function pna_uncertain(c)
    mag′ = abs(c)
    mag = mag′ ± mag_err(mag′)
    ang = rad2deg(angle(c)) ± phase_err(mag′)
    mag * cis(deg2rad(ang))
end

function read_to_df(fname)
    raw = DataFrame(CSV.File(fname))
    freqs = raw[!, :freqs] ./ 1e9
    amp_gamma = pna_uncertain.(raw[!, :amp_real] .+ raw[!, :amp_imag] .* 1im)
    cold_gamma = pna_uncertain.(raw[!, :cold_real] .+ raw[!, :cold_imag] .* 1im)
    hot_gamma = pna_uncertain.(raw[!, :hot_real] .+ raw[!, :hot_imag] .* 1im)
    sa_hot = col_meas(raw, "sa_hot")
    sa_cold = col_meas(raw, "sa_cold")
    dut_hot = col_meas(raw, "dut_hot")
    dut_cold = col_meas(raw, "dut_cold")
    tcold = tcold_meas.(freqs)
    DataFrame(
        freqs=freqs,
        hot_gamma=hot_gamma,
        amp_gamma=amp_gamma,
        cold_gamma=cold_gamma,
        dut_hot=dut_hot,
        dut_cold=dut_cold,
        sa_hot=sa_hot,
        sa_cold=sa_cold,
        tcold=tcold
    )
end

"""Compute the noise temperature of the device under test"""
function Tdut(ms, thot)
    # Compute Y factor
    hot = ms[!, :dut_hot] - ms[!, :sa_hot]
    cold = ms[!, :dut_cold] - ms[!, :sa_hot]
    Y = @. hot / cold

    # Compute correction
    Chot = correction.(ms[!, :hot_gamma], ms[!, :amp_gamma])
    Ccold = correction.(ms[!, :cold_gamma], ms[!, :amp_gamma])
    C = Chot ./ Ccold

    # Compute corrected noise temperature
    tcold = ms[!, :tcold]
    @. (thot - Y * tcold / C) / (Y / C - 1)
end

###### Data Analysis

function plot_summary(ms)
    f = Figure(; size=(800, 1000))
    ax = Axis(f[1, 1], title="Reflection Coefficients", xlabel="Freq (GHz)", ylabel="|S₁₁|² (dB)")
    lines!(ax, ms[!, :freqs] ./ 1e9, @. Measurements.value.(20 * log10(abs(ms[!, :hot_gamma]))); label="Hot Load")
    lines!(ax, ms[!, :freqs] ./ 1e9, @. Measurements.value.(20 * log10(abs(ms[!, :cold_gamma]))); label="Cold Load")
    lines!(ax, ms[!, :freqs] ./ 1e9, @. Measurements.value.(20 * log10(abs(ms[!, :amp_gamma]))); label="DUT")
    axislegend(ax, position=:rb)

    ax = Axis(f[2, 1], title="Powers", xlabel="Freq (GHz)", ylabel="dBm")
    lines!(ax, ms[!, :freqs] ./ 1e9, @. Measurements.value.(10 * log10(abs(ms[!, :hot]) / 1e-3)); label="Hot")
    lines!(ax, ms[!, :freqs] ./ 1e9, @. Measurements.value.(10 * log10(abs(ms[!, :cold]) / 1e-3)); label="Cold")
    lines!(ax, ms[!, :freqs] ./ 1e9, @. Measurements.value.(10 * log10(abs(ms[!, :sa]) / 1e-3)); label="RX")
    axislegend(ax, position=:rb)

    ax = Axis(f[3, 1], title="Y", xlabel="Freq (GHz)", ylabel="Y")
    hot = ms[!, :hot] - ms[!, :sa]
    cold = ms[!, :cold] - ms[!, :sa]
    Y = @. hot / cold
    lines!(ax, ms[!, :freqs] ./ 1e9, @. Measurements.value.(Y); label="Uncorrected")

    f
end

new_ms = read_to_df("$(@__DIR__)/../experimental_data/20240815_1444_ksWBLNAv2_SN1.csv")
old_ms = read_to_df("$(@__DIR__)/../experimental_data/20240815_1427_WBLNAv4_SN1.csv")
#plot_summary(old_ms)

# Thot using thermocouple error of 2C of 3sigma
Told = Tdut(old_ms, 297 ± 0.5)
Tnew = Tdut(new_ms, 297 ± 0.5)

function smooth_and_plot!(axis, T; label="")
    sg = savitzky_golay(Measurements.value.(T), 115, 1)
    sg_err = savitzky_golay(Measurements.uncertainty.(T), 115, 1)
    # Compute mean of unsmoothed data
    # Plot 2sigma
    mm = mean(T[50:700])
    μ = round(Measurements.value(mm), digits=2)
    σ = round(2 * Measurements.uncertainty(mm), digits=2)
    f = range(0.6, 2.1, 751)
    lines!(axis, f, sg.y, label="$label: $(μ) ± $(σ)")
    band!(axis, f, sg.y .- 2 * sg_err.y, sg.y .+ 2 * sg_err.y)
    # Save data
    pub_df = DataFrame(f=f, t=sg.y, e=2 .* sg_err.y)
    CSV.write("$(@__DIR__)/../publication_data/$(label).csv", pub_df)
    # Plot
    axis
end

begin
    f = Figure()
    ax = Axis(f[1, 1], title="Noise Temperature Comparison", xlabel="Freq (GHz)", ylabel="Noise (K)", limits=(nothing, (6, 16)))
    smooth_and_plot!(ax, Told; label="2-Step")
    smooth_and_plot!(ax, Tnew; label="NTL")
    axislegend(ax, position=:rb)
    f
end
