function updated_parameters!(m::MaxEnt)

    nJ = m.αJ * m.t^(-m.γJ)
    nh = m.αh * m.t^(-m.γh)

    Δx = nh * (m.x_obs .- m.x_mod)
    Δxy = nJ * (m.xy_obs .- m.xy_mod)

    m.h .+= Δx
    m.J .+= Δxy

    return nothing
end
