using MaxEntropyModel
using GLMakie
using Random, Distributions
using BenchmarkTools
using Statistics
using StatsBase
using LinearAlgebra

# %%  --------------------------------------------------------------------------
f = Figure()
ax = Axis(f[1, 1])
for q in 0.0:0.5:2
    lines!(ax, -5:0.01:5, x -> exp_q(-x^2, q), label="$q")
end
axislegend()
f
# %%  --------------------------------------------------------------------------


nspins = 16;
h = rand(Normal(0.0, 0.1), nspins);
J = rand(Normal(0.0, 0.1), nspins * (nspins - 1) ÷ 2);
# %%  --------------------------------------------------------------------------
S = map(x -> x < 0.5 ? -1 : 1, rand(1000, nspins));
m = MaxEnt(S, "teste", 'f');

m.h = copy(h);
m.J = copy(J);
full_tsallis!(m, 0.2);

E02 = copy(m.energy_hist);
full_tsallis!(m, 5.0);
E5 = copy(m.energy_hist);

# %%  --------------------------------------------------------------------------
nbins = 64
h02 = fit(Histogram, E02, nbins=nbins)
h02 = normalize(h02, mode=:pdf)
x02 = h02.edges[1]
x02 = 0.5 * (x02[2:end] + x02[1:end-1])

h5 = fit(Histogram, E5, nbins=nbins)
h5 = normalize(h5, mode=:pdf)
x5 = h5.edges[1]
x5 = 0.5 * (x5[2:end] + x5[1:end-1])


f = Figure()
ax = Axis(f[1, 1])
stairs!(ax, x02, h02.weights)
stairs!(ax, x5, h5.weights)
f
