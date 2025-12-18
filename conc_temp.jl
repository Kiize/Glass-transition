using DelimitedFiles

data = readdlm("mean_concentration.txt", skipstart=1)
c_mean, β, L = data[:, 1], data[:, 2], data[:, 3]
scatter(β, c_mean, plot_title="Concentration vs temperature", xlabel=L"$\beta$", ylabel=L"C", legend=false)
savefig("glass_concentration_vs_temperature.png")