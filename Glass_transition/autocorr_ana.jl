using Plots
using ProgressMeter
using LaTeXStrings
using DelimitedFiles
using LsqFit

β_arr = [0.001, 0.01, 0.05, 0.1, 0.25, 0.4, 0.5, 0.6, 0.75, 0.8, 1.0, 1.25, 1.5,
        2.0, 5.0, 10.0]

# Fit

model(t, p) = p[1] .* exp.(-(t./p[2])) # p[1] is amplitude, p[2] is tau
p0 = [1.0, 2.0]  # Initial guess for amplitude, tau

open("autocorr_rel_time.txt", "w") do io
    writedlm(io, ["tau" "beta"])
end

@showprogress for β in β_arr
    local f = "autocorr_beta$(replace(string(β), "." => "-")).txt"
    local y = readdlm(f, skipstart=1)
    local t = range(start = 0, length = length(y))

    # Perform the curve fit
    local fit = curve_fit(model, t, y[:], p0)

    # Extract fitted parameters
    local fitted_params = fit.param
    local tau = fitted_params[2]

    open("autocorr_rel_time.txt", "a") do io
        writedlm(io, [tau β])
    end
end

# Relaxation time vs temperature

data = readdlm("autocorr_rel_time.txt", skipstart=1)
τ, β_arr = data[:, 1], data[:, 2]
scatter(β_arr, τ, legend=false, plot_title="Relaxation time vs temperature", ylabel=L"$\tau$", xlabel=L"$\beta$")
savefig("rel_time_vs_temp.png")
scatter(β_arr, τ, legend=false, plot_title="Relaxation time vs temperature in semilog scale", ylabel=L"$\tau$", xlabel=L"$\beta$", yscale=:log10)
savefig("rel_time_vs_temp_log.png")
