using Base: nothing_sentinel
using MaxEntropyModel
using GLMakie
using Random, Distributions
using BenchmarkTools
using Statistics
using StatsBase
using LinearAlgebra

# %%
# generate model params
nspins = 16;
h = rand(Normal(0.0, 0.1), nspins);
h = 2 * rand(nspins) .- 1
J = rand(Normal(0.0, 0.2), nspins * (nspins - 1) ÷ 2);

# %%

reg = false

S = map(x -> x < 0.5 ? -1 : 1, rand(1000, nspins));
m1 = MaxEnt(S, "teste", 'f');


m1.h = copy(h);
m1.J = copy(J);
m1.q = 1.0

full_measurements!(m1);

# %%  --------------------------------------------------------------------------
m2 = MaxEnt(m1, "teste2", 'q')
m2.q = 0.9
m2.tol = 1.0e-7
m2.n_relax_steps = 2000
max_entropy_relax!(m2)

# %%

nbins = 150
f = Figure()
ax1 = Axis(f[1, 1], xlabel="E", ylabel="N(E)")
ax2 = Axis(f[1, 1], ylabel="P(E)", yaxisposition=:right, rightspinevisible=false)

hidexdecorations!(ax2, grid=true)
hidespines!(ax2)
linkxaxes!(ax1, ax2)

m = m2;
q = m.q;
H0 = m.H0_vals[1];
var = reg ? m.Hj_vals .+ H0 : m.Hj_vals;
vmax = maximum(var);
vmin = minimum(var);
println((vmax, vmin))
dx = (vmax - vmin) / nbins;
edges = LinRange(vmin - dx / 2, vmax + dx / 2, nbins + 1);
hvar = fit(Histogram, var, edges);
y = hvar.weights;
x = hvar.edges[1];
x = 0.5 * (x[2:end] + x[1:end-1]);

Pj = zeros(size(x))
for j in eachindex(var)
    hj = var[j]
    k = findfirst(x -> x > hj, edges) - 1
    Pj[k] += m.Pj_vals[j]
end
Pj /= sum(Pj)
println(Pj)

P = reg ? y .* exp_q.(-(x + H0, q)) : y .* exp_q.(-x, q);
Z = sum(P);
y1 = exp_q.(-(x), q);
y1 ./= sum(y1);
y1 /= dx;
lines!(ax1, x, y, label="N(E)", color=:dodgerblue);
# lines!(ax2, x, P / dx / Z, linewidth=2, color=:black, label="P(E)");
lines!(ax2, x, Pj / dx, linewidth=2, color=:black, label="P(E)");
lines!(ax2, x, y1, color=:red, linestyle=:dash, label="exp_q(-E)");
Label(f[1, 1, Top()], "q = $q", fontsize=24, padding=(0, 0, 10, 0));
axislegend(ax1, position=:lt, framevisible=false);
axislegend(ax2, position=:rt, framevisible=false);
f

# %%








# %%
# just test exp_q
f = Figure();
ax1 = Axis(f[1, 1]);
for q in 0.0:0.5:2
    lines!(ax1, -5:0.01:5, x -> exp_q(-x^2, q), label="$q")
end
axislegend();
f

