function energy(m::MaxEnt)
    H = -dot(m.h, m.s)

    t = 1
    for i in 1:m.nspins-1, j in i+1:m.nspins
        @inbounds H -= m.s[i] * m.J[t] * m.s[j]
        t += 1
    end

    return H
end


function deltaEnergy(m::MaxEnt, i::Int64)
    ΔH = 2 * m.s[i] * m.h[i]

    for j in 1:m.nspins
        if j != i
            t = m.bond[i, j]
            ΔH += 2.0 * m.s[i] * m.J[t] * m.s[j]
        end
    end

    return ΔH
end
