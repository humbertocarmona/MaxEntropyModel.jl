
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

        m.x_mod = m.x_mod + Zx * m.s
        t = 1
        for i in 1:nspins-1, j in i+1:nspins
            m.xy_mod[t] += m.s[i] * m.s[j] * Zx
            t += 1
        end
    end
    m.x_mod = m.x_mod / Z
    m.xy_mod = m.xy_mod / Z

    return nothing
end


