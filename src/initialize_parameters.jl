function init_parameters!(m::MaxEnt, h::Vector{Float64}, J::Vector{Float64})
    @assert size(h) == size(m.x_obs) "expected size(h)=$(m.x_obs), found $(size(h))"
    @assert size(J) == size(m.xy_mod) "expected lengh(J)=$(m.xy_mod), found $(size(J))"

    m.h = h
    m.J = J

end

function init_parameters!(m::MaxEnt, μh::Float64, σh::Float64, μJ::Float64, σJ::Float64)
    disth = Normal(μh, σh)
    distJ = Normal(μJ, σJ)
    m.h = rand(disth, size(m.x_obs))
    m.J = rand(distJ, size(m.xy_obs))
end

function init_parameters!(m::MaxEnt)
    m.h = copy(m.x_obs)
    m.J = copy(m.xy_obs)
end
