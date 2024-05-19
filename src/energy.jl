function energy(m::MaxEnt)

    H = -dot(m.h, m.sj)

    t = 1
    for j in 1:m.nspins-1
        for i in j+1:m.nspins
            @inbounds H -= m.sj[i] * m.J[t] * m.sj[j]
            t += 1
        end
    end

    return H
end


# function deltaEnergy(m::MaxEnt, i::Int64)
#     ΔH = 2 * m.sj[i] * m.h[i]

#     for j in 1:m.nspins
#         if j != i
#             t = m.bond[i, j]
#             ΔH += 2.0 * m.sj[i] * m.J[t] * m.sj[j]
#         end
#     end

#     return ΔH
# end

# should always prefer to loop vertically in Julia
function deltaEnergy(m::MaxEnt, j::Int64)
    ΔH = 2 * m.sj[j] * m.h[j]

    for i in 1:m.nspins
        if i != j
            t = m.bond[i, j]
            ΔH += 2.0 * m.sj[j] * m.J[t] * m.sj[j]
        end
    end

    return ΔH
end