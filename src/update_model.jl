function updated_parameters!(m::MaxEnt)

    nJ = m.ηJ * m.t^(-m.γJ)
    nh = m.ηh * m.t^(-m.γh)

    Δx = nh * (m.x_obs .- m.x_mod)
    Δxy = nJ * (m.xy_obs .- m.xy_mod)

    m.h .+= Δx + m.α * m.Δx
    m.J .+= Δxy + m.α * m.Δxy

    m.Δxy = copy(Δxy)
    m.Δx = copy(Δx)

    return nothing
end


