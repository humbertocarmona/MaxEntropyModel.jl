function full_iteration!(model::MaxEnt)
	nspins = model.nspins
	@assert nspins < 25 "maximum full ensemble system is 24, found $nspins"

	Z = 0.0
	model.x_mod .= zeros(nspins)
	model.xy_mod .= zeros(Float64, nspins * (nspins - 1) ÷ 2)
	j = 1
	s = zeros(Int, nspins)
	ps = zeros(Float64, nspins)

	for p in spin_permutations_iterator(nspins)
		s .= collect(p)
		Hj = energy(model, s)
		Pj = exp(-model.β * Hj)
		ps .= Pj * s
		Z += Pj

		model.x_mod .= model.x_mod .+ ps
		i = 1
		for k in 1:nspins-1
			for l in k+1:nspins
				@inbounds model.xy_mod[i] += s[k] * ps[l]
				i += 1
			end
		end
		model.Hj_vals[j] = Hj
		model.Pj_vals[j] = Pj

		j += 1
	end

	model.x_mod ./= Z
	model.xy_mod ./= Z
	model.Pj_vals ./= Z
	return nothing
end

function full_measurements!(model::MaxEnt)
	nspins = model.nspins
	@assert nspins < 25 "maximum full ensemble system is 24, found $nspins"

	Z = 0.0
	model.x_mod .= zeros(nspins)
	model.xy_mod .= zeros(Float64, nspins * (nspins - 1) ÷ 2)
	model.Hj_vals = Array{Float64}(undef, 2^model.nspins)

	model.H_mean = 0.0
	H2_mean = 0.0
	model.M_mean = 0.0
	model.xyz_mod .= zeros(nspins * (nspins - 1) * (nspins - 2) ÷ 6)
	model.ones_dist_mod .= zeros(nspins + 1)

	j = 1
	s = zeros(Int, nspins)
	ps = zeros(Float64, nspins)

	for p in spin_permutations_iterator(nspins)
		s .= collect(p)
		Hj = energy(model, s)
		Pj = exp(-model.β * Hj)
		ps .= Pj * s

		Z += Pj

		model.x_mod .= model.x_mod .+ ps
		i = 1
		for k in 1:nspins-1, l in k+1:nspins
			model.xy_mod[i] += s[k] * ps[l]
			i += 1
		end

		model.H_mean += Hj * Pj
		H2_mean += Hj * Hj * Pj
		model.M_mean += sum(s) * Pj
		i = 1
		for k in 1:nspins-2
			for l in k+1:nspins-1
				for m in l+1:nspins
					model.xyz_mod[i] += s[k] * s[l] * ps[m]
					i += 1
				end
			end
		end
		k = count(isone.(s))
		model.ones_dist_mod[k+1] += Pj

		model.Hj_vals[j] = Hj
		model.Pj_vals[j] = Pj
		j += 1
	end
	model.Pj_vals ./= Z
	model.x_mod ./= Z
	model.xy_mod ./= Z
	model.H_mean /= Z
	H2_mean /= Z
	model.CV = model.β^2 * (H2_mean - model.H_mean^2)
	model.M_mean /= Z
	model.xyz_mod ./= Z
	model.ones_dist_mod ./= sum(model.ones_dist_mod)

	# pearson_mod!(model)
	return nothing
end
