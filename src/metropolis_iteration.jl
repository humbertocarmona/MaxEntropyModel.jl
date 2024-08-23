function flip!(s::Vector{Int}, i::Int64, H::Float64, model::MaxEnt, rng)
    ΔE = deltaEnergy(model, s, i)
    if ΔE < 0.0 || rand(rng) < exp(-model.β * ΔE)
        s[i] = -s[i]
        H += ΔE
    end
    return H
end

function metropolis_iteration!(model::MaxEnt)
    rng = Xoshiro(model.mc_seed)
    nspins = model.nspins

    s = map(x -> x < 0.5 ? -1 : 1, rand(nspins))
    Hj = energy(model, s)

    model.x_mod .= zeros(nspins)
    model.xy_mod .= zeros(Float64, nspins * (nspins - 1) ÷ 2)
    model.Pj_vals = Array{Float64}(undef, model.n_samples)
    model.Hj_vals = Array{Float64}(undef, model.n_samples)

    n_flips = model.n_samples * model.n_coherence / model.n_rept + model.n_equilibrium

    j = 1
    for _ in 1:model.n_rept
        for n in 1:n_flips
            t = n - model.n_equilibrium

            i = rand(rng, 1:model.nspins)
            Hj = flip!(s, i, Hj, model, rng)

            if (t > 0) && (t % model.n_coherence == 0)
                model.x_mod .= model.x_mod .+ s
                i = 1
                for k in 1:nspins-1, l in k+1:nspins
                    model.xy_mod[i] += s[k] * s[l]
                    i += 1
                end
                model.Hj_vals[j] = Hj
                j += 1
            end
        end
    end
    # s -= 1
    # @assert s == m.n_samples "expected s=$(m.n_samples), got $s"

    model.x_mod ./= model.n_samples
    model.xy_mod ./= model.n_samples

    return nothing
end

function metropolis_measurements!(model::MaxEnt)
    rng = Xoshiro(model.mc_seed)
    nspins = model.nspins

    s = map(x -> x < 0.5 ? -1 : 1, rand(nspins))
    Hj = energy(model, s)

    model.x_mod .= zeros(nspins)
    model.xy_mod .= zeros(Float64, nspins * (nspins - 1) ÷ 2)
    model.Pj_vals = Array{Float64}(undef, model.n_samples)
    model.Hj_vals = Array{Float64}(undef, model.n_samples)

    model.H_mean = 0.0
    H2_mean = 0.0
    model.M_mean = 0.0
    model.xyz_mod .= zeros(nspins * (nspins - 1) * (nspins - 2) ÷ 6)
    model.ones_dist_mod .= zeros(nspins + 1)
    samples = zeros(Int64, (model.n_samples, nspins))

    n_flips = model.n_samples * model.n_coherence / model.n_rept + model.n_equilibrium
    j = 1
    for _ in 1:model.n_rept
        for n in 1:n_flips
            t = n - model.n_equilibrium

            i = rand(rng, 1:model.nspins)
            Hj = flip!(s, i, Hj, model, rng)

            if (t > 0) && (t % model.n_coherence == 0)
                model.x_mod .= model.x_mod .+ s
                i = 1
                for k in 1:nspins-1, l in k+1:nspins
                    model.xy_mod[i] += s[k] * s[l]
                    i += 1
                end
                model.Hj_vals[j] = Hj
                model.H_mean += Hj
                H2_mean += Hj * Hj
                model.M_mean += sum(s)
                i = 1
                for k in 1:nspins-2
                    for l in k+1:nspins-1
                        for m in l+1:nspins
                            model.xyz_mod[i] += s[k] * s[l] * s[m]
                            i += 1
                        end
                    end
                end
                k = count(isone.(s))
                model.ones_dist_mod[k+1] += 1
                samples[j, :] = copy(s)
                j += 1
            end
        end
    end
    # s -= 1
    # @assert s == m.n_samples "expected s=$(m.n_samples), got $s"

    model.x_mod ./= model.n_samples
    model.xy_mod ./= model.n_samples
    # pearson_mod!(model)
    model.H_mean /= model.n_samples
    H2_mean /= model.n_samples
    model.CV = model.β^2 * (H2_mean - model.H_mean^2)
    model.M_mean /= model.n_samples
    model.xyz_mod ./= model.n_samples
    model.ones_dist_mod ./= sum(model.ones_dist_mod)
    return samples
end

function metropolis_measurements!(model::MaxEnt, samples::Matrix{Int64})
    nsamples, nspins = size(samples)
    @assert model.nspins == nspins "got m.nspins=$(model.nspins) != $(nspins)"

    model.x_mod .= zeros(nspins)
    model.xy_mod .= zeros(Float64, nspins * (nspins - 1) ÷ 2)
    model.Hj_vals = Array{Float64}(undef, model.n_samples)
    model.H_mean = 0.0
    H2_mean = 0.0
    model.M_mean = 0.0
    model.xyz_mod .= zeros(nspins * (nspins - 1) * (nspins - 2) ÷ 6)
    model.ones_dist_mod .= zeros(nspins + 1)

    for j in 1:nsamples
        s = samples[:, j]
        Hj = energy(model, s)
        model.x_mod .= model.x_mod .+ s
        i = 1
        for k in 1:nspins-1, l in k+1:nspins
            model.xy_mod[i] += s[k] * s[l]
            i += 1
        end
        model.H_mean += Hj
        H2_mean += Hj * Hj
        model.M_mean += sum(sj)
        i = 1
        for k in 1:nspins-2
            for l in k+1:nspins-1
                for m in l+1:nspins
                    model.xyz_mod[i] += s[k] * s[l] * s[m]
                    i += 1
                end
            end
        end
        k = count(isone.(s))
        model.ones_dist_mod[k+1] += 1
        model.Hj_vals[j] = Hj
    end

    model.x_mod ./= model.n_samples
    model.xy_mod ./= model.n_samples
    model.H_mean /= model.n_samples
    H2_mean /= model.n_samples
    model.CV = model.β^2 * (H2_mean - model.H_mean^2)
    model.M_mean /= model.n_samples
    model.xyz_mod ./= model.n_samples
    model.ones_dist_mod ./= sum(model.ones_dist_mod)

    # pearson_mod!(model)
end
