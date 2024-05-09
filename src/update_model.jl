function updated_parameters!(m::MaxEnt, pearson=true)

    nJ = m.α * m.t^(-m.γ)
    nh = 2.0 * nJ

    Δx = m.x_obs .- m.x_mod
    if pearson
        Δxy = m.pearson_obs .- m.pearson_mod
    else
        Δxy = m.xy_obs .- m.xy_mod
    end

    m.h .= m.h .+ nh * Δx
    m.J .= m.J .+ nJ * Δxy

    return nothing
end
