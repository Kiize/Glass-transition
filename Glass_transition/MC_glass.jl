using Random
using Plots
using ProgressMeter
using LaTeXStrings
using StatsBase
using DelimitedFiles

#= 
Initialize the LxL lattice with random variables which can take +-1 values. 
=#
function init_lattice(L)
    lattice = rand([0, 1], L, L)
    return lattice
end

function delta_energy(lattice, i, j)
    return -2 * lattice[i, j] + 1
end

function check_neigh(lattice, i, j)
    L = size(lattice, 1)

    # Check with PBC

    return lattice[mod1(i+1, L), j] == 1 ||
           lattice[mod1(i-1, L), j] == 1 ||
           lattice[i, mod1(j+1, L)] == 1 ||
           lattice[i, mod1(j-1, L)] == 1

end

function metropolis_step!(lattice, β, n_therm, heatmap)
    L = size(lattice, 1)

    for _ in 1:n_therm
        # Choose a random site (i, j)
        i = rand(1:L)
        j = rand(1:L)

        # Check if at least one neighbour is 1, skip otherwise
        !check_neigh(lattice, i, j) && continue

        ΔE = delta_energy(lattice, i, j)

        if ΔE <= 0 || rand() < exp(-β * ΔE)
            lattice[i, j] = mod(lattice[i, j] + 1, 2)
            heatmap[i, j] += 1
        end
    end

    return sum(lattice)/L^2
end

# Simulazione

L = 128
β_arr = [0.5]
T = 10_000
MC_time = L^2
n_therm = 5000

concentration = zeros(length(β_arr), T)
lattice = init_lattice(L)

# Uncomment to re-write the entire file
#= open("mean_concentration.txt", "w") do io
    writedlm(io, ["mean" "beta" "L"])
end =#

for (i, β) in enumerate(β_arr)
    heat_map = zeros(L, L)

    # Thermalization

    @showprogress for t in 1:n_therm
        _ = metropolis_step!(lattice, β, MC_time, heat_map)
    end

    # Measuring concentration
    
    heat_map = zeros(L, L) # heatmap reset

    @showprogress for t in 1:T
        concentration[i, t] = metropolis_step!(lattice, β, MC_time, heat_map)
    end

    #= open("mean_concentration.txt", "a") do io
        writedlm(io, [mean(concentration[i, :]) β L])
    end =#

    c_autocor = autocor(concentration[i, :])

    #= open("autocorr_beta$(replace(string(β), "." => "-")).txt", "w") do io
        writedlm(io, c_autocor)
    end =#

    heatmap(heat_map, color=:RdBu, plot_title = L"Heatmap at $L$ = %$(L) and $\beta$ = %$(β)") 
    savefig("glass_heatmap_beta$(replace(string(β), "." => "-")).png")
end


#= heatmap(lattice, color=:RdBu, clims=(0,1), plot_title = L"Lattice at $L$ = %$(L) and $\beta$ = %$(β)") # 
savefig("glass_lattice_L$(L)_beta$(replace(string(β), "." => "-")).png")
scatter(concentration, plot_title = L"Concentration at $L$ = %$(L) and $\beta$ = %$(β)")
savefig("glass_lattice_concentration_L$(L)_beta$(replace(string(β), "." => "-")).png")
c_autocor = autocor(concentration)
scatter(c_autocor, plot_title = L"Autocorrelation at $L$ = %$(L) and $\beta$ = %$(β)")
savefig("glass_lattice_autocorr_L$(L)_beta$(replace(string(β), "." => "-")).png") =#


