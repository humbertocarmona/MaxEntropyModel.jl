function energy(m::MaxEnt, s::Vector{Int})
	H = -dot(m.h, s)

	t = 1
	@inbounds for j in 1:m.nspins-1
		for i in j+1:m.nspins
			H -= s[i] * m.J[t] * s[j]
			t += 1
		end
	end

	return H
end

# should always prefer to loop vertically in Julia
# TODO: checar essa soma!!!!
function deltaEnergy(m::MaxEnt, s::Vector{Int}, j::Int64)
	ΔH = 2 * s[j] * m.h[j]

	for i in 1:m.nspins
		if i != j
			t = m.bond[i, j]
			ΔH += 2.0 * s[i] * m.J[t] * s[j] 
		end
	end

	return ΔH
end
