using Base: nothing_sentinel
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


nspins = 20;
h = rand(Normal(0.0, 0.1), nspins);
# h = 2 * rand(nspins) .- 1
J = rand(Normal(0.0, 0.1), nspins * (nspins - 1) ÷ 2);
# %%  --------------------------------------------------------------------------
S = map(x -> x < 0.5 ? -1 : 1, rand(1000, nspins));
m1 = MaxEnt(S, "teste", 'f');
m2 = MaxEnt(S, "teste", 'f');

m1.h = copy(h);
m1.J = copy(J);
m2.h = copy(h);
m2.J = copy(J);

q1 = 5.0
full_tsallis!(m1, q1);

q2 = 0.0
full_tsallis!(m2, q2);
# %%  --------------------------------------------------------------------------
f = Figure()
ax = Axis(f[1, 1], yscale=log10)
ylims!(ax, 1e-9, 1)

x1 = (m1.PE_edges[1:end-1] + m1.PE_edges[2:end]) / 2.0
x2 = (m2.PE_edges[1:end-1] + m2.PE_edges[2:end]) / 2.0
y1 = 1.0 * m1.PE_weights
y2 = 1.0 * m2.PE_weights

lines!(ax, x1, y1, label="$(q1)")
lines!(ax, x2, y2, label="$(q2)")
axislegend()
f


# %%  --------------------------------------------------------------------------
nbins = 204
E1 = 1 .- (1 - q1) .* (m1.H_vals .+ m1.H0_vals[1])
h1 = fit(Histogram, E1, nbins=nbins)
h1 = normalize(h1, mode=:pdf)
x1 = h1.edges[1]
x1 = 0.5 * (x1[2:end] + x1[1:end-1])

E2 = 1 .- (1 - q2) .* (m2.H_vals .+ m2.H0_vals[1])
h2 = fit(Histogram, E2, nbins=nbins)
h2 = normalize(h2, mode=:pdf)
x2 = h2.edges[1]
x2 = 0.5 * (x2[2:end] + x2[1:end-1])


f = Figure()
ax = Axis(f[1, 1])
xlims!(ax, 0, nothing)
scatter!(ax, x1, h1.weights)
scatter!(ax, x2, h2.weights)
f


