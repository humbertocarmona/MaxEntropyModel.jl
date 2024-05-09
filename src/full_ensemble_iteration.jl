
function full_iteration!(m::MaxEnt)
    nspins = m.nspins
    @assert nspins < 25 "maximum full enesemble system is 24, found $nspins"

    Z = 0.0
    m.x_mod .= zeros(nspins)

    m.xy_mod .= zeros(Float64, nspins * (nspins - 1) ÷ 2)


    for p in spin_permutations_iterator(nspins)
        m.s .= collect(p)
        Ex = energy(m)
        Zx = exp(-m.β * Ex)
        Z += Zx
        m.x_mod .= m.x_mod .+ Zx .* m.s
        t = 1
        for i in 1:nspins-1, j in i+1:nspins
            m.xy_mod[t] += m.s[i] * m.s[j] * Zx
            t += 1
        end
    end
    m.x_mod ./= Z
    m.xy_mod ./= Z

    pearson_mod!(m)

    return nothing
end


function full_iteration!(m::MaxEnt, meas::Bool)
    nspins = m.nspins
    @assert nspins < 25 "maximum full enesemble system is 24, found $nspins"

    Z = 0.0
    m.x_mod .= zeros(nspins)
    m.xy_mod .= zeros(Float64, nspins * (nspins - 1) ÷ 2)
    if meas
        m.energy_mean = 0.0
        m.specific_heat = 0.0
        m.magnetization_mean = 0.0
        m.xyz_mod .= zeros(nspins * (nspins - 1) * (nspins - 2) ÷ 6)
        m.ones_dist_mod .= zeros(nspins + 1)
    end

    for p in spin_permutations_iterator(nspins)
        m.s .= collect(p)
        Ex = energy(m)
        Zx = exp(-m.β * Ex)
        Z += Zx

        m.x_mod .= m.x_mod .+ Zx .* m.s
        t = 1
        for i in 1:nspins-1, j in i+1:nspins
            m.xy_mod[t] += m.s[i] * m.s[j] * Zx
            t += 1
        end

        if meas
            m.energy_mean += Ex * Zx
            m.specific_heat += Ex * Ex * Zx
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
        end

    end
    m.x_mod ./= Z
    m.xy_mod ./= Z

    if meas
        m.energy_mean /= Z
        m.specific_heat /= Z
        m.specific_heat = m.β^2 * (m.specific_heat - m.energy_mean^2)
        m.magnetization_mean /= Z
        m.xyz_mod ./= Z
        m.ones_dist_mod ./= sum(m.ones_dist_mod)
    end

    pearson_mod!(m)
    return nothing
end