function full_tsallis!(m::MaxEnt, q::Float64)
    nspins = m.nspins
    @assert nspins < 25 "maximum full enesemble system is 24, found $nspins"

    debug(LOGGER, "aqui q = $q")
    Zq = 0.0
    m.x_mod .= zeros(nspins)
    m.xy_mod .= zeros(Float64, nspins * (nspins - 1) ÷ 2)

    H0 = compute_energy_shift(m, q)
    m.H_vals = Array{Float64}(undef, 2^m.nspins)
    m.H0_vals[m.t] = H0

    m.PE_weights .= 0.0

    s = 1
    for p in spin_permutations_iterator(nspins)
        m.s .= collect(p)
        m.H = energy(m)
        Zx = exp_q(-m.β * (m.H + H0), q)
        Zq += Zx
        m.H_vals[s] = m.H
        m.x_mod .= m.x_mod .+ Zx .* m.s
        t = 1
        for i in 1:nspins-1
            for j in i+1:nspins
                @inbounds m.xy_mod[t] += m.s[i] * m.s[j] * Zx
                t += 1
            end
        end
        if m.Hmin <= m.H <= m.Hmax
            k = findfirst(x -> x > m.H, m.PE_edges) - 1
            m.PE_weights[k] += Zx
        end
        s += 1
    end

    m.PE_weights ./= sum(m.PE_weights)

    m.x_mod ./= Zq
    m.xy_mod ./= Zq

    return nothing
end

function full_tsallis_measurements!(m::MaxEnt, q::Float64)
    nspins = m.nspins
    @assert nspins < 25 "maximum full enesemble system is 24, found $nspins"

    Zq = 0.0
    m.x_mod .= zeros(nspins)
    m.xy_mod .= zeros(Float64, nspins * (nspins - 1) ÷ 2)
    m.H_vals = Array{Float64}(undef, 2^m.nspins)

    m.H_mean = 0.0
    H2_mean = 0.0
    m.magnetization_mean = 0.0
    m.xyz_mod .= zeros(nspins * (nspins - 1) * (nspins - 2) ÷ 6)
    m.ones_dist_mod .= zeros(nspins + 1)

    H0 = compute_energy_shift(m, q)
    m.H0_vals[m.t] = H0

    s = 1
    for p in spin_permutations_iterator(nspins)
        m.s .= collect(p)
        m.H = energy(m)
        Zx = exp_q(-m.β * (m.H + H0), q)
        Zq += Zx

        m.x_mod .= m.x_mod .+ Zx .* m.s
        t = 1
        for i in 1:nspins-1, j in i+1:nspins
            m.xy_mod[t] += m.s[i] * m.s[j] * Zx
            t += 1
        end

        m.H_mean += m.H * Zx
        H2_mean += m.H * m.H * Zx
        m.magnetization_mean += sum(m.s) * Zx
        t = 1
        for i in 1:nspins-2
            for j in i+1:nspins-1
                for k in j+1:nspins
                    m.xyz_mod[t] += m.s[i] * m.s[j] * m.s[k] * Zx
                    t += 1
                end
            end
        end
        k = count(isone.(m.s))
        m.ones_dist_mod[k+1] += Zx

        m.H_vals[s] = m.H
        s += 1
    end
    m.x_mod ./= Zq
    m.xy_mod ./= Zq
    m.H_mean /= Zq
    H2_mean /= Zq
    m.specific_heat = m.β^2 * (H2_mean - m.H_mean^2)
    m.magnetization_mean /= Zq
    m.xyz_mod ./= Zq
    m.ones_dist_mod ./= sum(m.ones_dist_mod)

    pearson_mod!(m)
    return nothing
end

