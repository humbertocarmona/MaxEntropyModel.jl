function flip!(m::MaxEnt, i::Int64, rng)
    ΔE = deltaEnergy(m, i)
    if ΔE < 0.0 || rand(rng) < exp(-m.β * ΔE)
        m.s[i] = -m.s[i]
        m.Es += ΔE
    end
end

function metropolis_iteration!(m::MaxEnt)
    rng = Xoshiro(m.mc_seed)
    nspins = m.nspins

    m.Es = energy(m)

    m.x_mod .= zeros(nspins)
    m.xy_mod .= zeros(Float64, nspins * (nspins - 1) ÷ 2)

    n_flips = m.n_samples * m.n_coherence / m.n_rept + m.n_equilibrium

    s = 1
    for _ in 1:m.n_rept
        m.s .= -1
        for n in 1:n_flips
            t = n - m.n_equilibrium

            s_flip = rand(rng, 1:m.nspins)
            flip!(m, s_flip, rng)

            if (t > 0) && (t % m.n_coherence == 0)
                m.x_mod .= m.x_mod .+ m.s
                k = 1
                for i in 1:nspins-1, j in i+1:nspins
                    m.xy_mod[k] += m.s[i] * m.s[j]
                    k += 1
                end
                s += 1
            end
        end
    end
    # s -= 1
    # @assert s == m.n_samples "expected s=$(m.n_samples), got $s"

    m.x_mod ./= m.n_samples
    m.xy_mod ./= m.n_samples

    return nothing
end


function metropolis_measurements!(m::MaxEnt)
    rng = Xoshiro(m.mc_seed)
    nspins = m.nspins

    m.Es = energy(m)

    m.x_mod .= zeros(nspins)
    m.xy_mod .= zeros(Float64, nspins * (nspins - 1) ÷ 2)
    m.energy_hist = Array{Float64}(undef, m.n_samples)
    m.energy_mean = 0.0
    m.specific_heat = 0.0
    m.magnetization_mean = 0.0
    m.xyz_mod .= zeros(nspins * (nspins - 1) * (nspins - 2) ÷ 6)
    m.ones_dist_mod .= zeros(nspins + 1)
    samples = zeros(Int64, (m.n_samples, nspins))

    n_flips = m.n_samples * m.n_coherence / m.n_rept + m.n_equilibrium

    s = 1
    for _ in 1:m.n_rept
        m.s .= -1
        for n in 1:n_flips
            t = n - m.n_equilibrium

            s_flip = rand(rng, 1:m.nspins)
            flip!(m, s_flip, rng)

            if (t > 0) && (t % m.n_coherence == 0)
                m.x_mod .= m.x_mod .+ m.s
                k = 1
                for i in 1:nspins-1, j in i+1:nspins
                    m.xy_mod[k] += m.s[i] * m.s[j]
                    k += 1
                end
                m.energy_hist[s] = m.Es
                m.energy_mean += m.Es
                m.specific_heat += m.Es * m.Es
                m.magnetization_mean += sum(m.s)
                t = 1
                for i in 1:nspins-2
                    for j in i+1:nspins-1
                        for k in j+1:nspins
                            m.xyz_mod[t] += m.s[i] * m.s[j] * m.s[k]
                            t += 1
                        end
                    end
                end
                k = count(isone.(m.s))
                m.ones_dist_mod[k+1] += 1
                samples[s, :] = copy(m.s)
                s += 1
            end
        end
    end
    # s -= 1
    # @assert s == m.n_samples "expected s=$(m.n_samples), got $s"

    m.x_mod ./= m.n_samples
    m.xy_mod ./= m.n_samples
    pearson_mod!(m)
    m.energy_mean /= m.n_samples
    m.specific_heat /= m.n_samples
    m.specific_heat = m.β^2 * (m.specific_heat - m.energy_mean^2)
    m.magnetization_mean /= m.n_samples
    m.xyz_mod ./= m.n_samples
    m.ones_dist_mod ./= sum(m.ones_dist_mod)
    return samples
end


function metropols_measurements!(m::MaxEnt, samples::Matrix{Int64})
    nsamples, nspins = size(samples)
    @assert m.nspins == nspins "got m.nspins=$(m.nspins) != $(nspins)"

    m.x_mod .= zeros(nspins)
    m.xy_mod .= zeros(Float64, nspins * (nspins - 1) ÷ 2)
    m.energy_hist = Array{Float64}(undef, m.n_samples)
    m.energy_mean = 0.0
    m.specific_heat = 0.0
    m.magnetization_mean = 0.0
    m.xyz_mod .= zeros(nspins * (nspins - 1) * (nspins - 2) ÷ 6)
    m.ones_dist_mod .= zeros(nspins + 1)

    for n in 1:nsamples
        m.s = samples[:, n]
        m.x_mod .= m.x_mod .+ m.s
        k = 1
        for i in 1:nspins-1, j in i+1:nspins
            m.xy_mod[k] += m.s[i] * m.s[j]
            k += 1
        end
        m.energy_mean += m.Es
        m.specific_heat += m.Es * m.Es
        m.magnetization_mean += sum(m.s)
        t = 1
        for i in 1:nspins-2
            for j in i+1:nspins-1
                for k in j+1:nspins
                    m.xyz_mod[t] += m.s[i] * m.s[j] * m.s[k]
                    t += 1
                end
            end
        end
        k = count(isone.(m.s))
        m.ones_dist_mod[k+1] += 1
        m.energy_hist[s] = m.Es
    end

    m.x_mod ./= m.n_samples
    m.xy_mod ./= m.n_samples
    m.energy_mean /= m.n_samples
    m.specific_heat /= m.n_samples
    m.specific_heat = m.β^2 * (m.specific_heat - m.energy_mean^2)
    m.magnetization_mean /= m.n_samples
    m.xyz_mod ./= m.n_samples
    m.ones_dist_mod ./= sum(m.ones_dist_mod)

    pearson_mod!(m)


end
