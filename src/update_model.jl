function updated_parameters!(m::MaxEnt)

    nJ = m.α * m.t^(-m.γ)
    nh = 2.0 * nJ

    Δx = nh * (m.x_obs .- m.x_mod)
    Δxy = nJ * (m.xy_obs .- m.xy_mod)

    m.h .= m.h .+ Δx
    m.J .= m.J .+ Δxy

    return nothing
end
