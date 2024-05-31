
"""
	full_q_iteration!(model::MaxEnt)

	this uses exp_q, that returns exp if q==1
"""
function full_q_iteration!(model::MaxEnt)
	nspins = model.nspins
	@assert nspins < 25 "maximum full ensemble system is 24, found $nspins"

	q = model.q
	Zq = 0.0

	model.x_mod .= zeros(nspins)
	model.xy_mod .= zeros(Float64, nspins * (nspins - 1) ÷ 2)
	
	H0 = 0.0
	if model.reg
		H0 = compute_energy_shift(model, q)
		model.H0_vals[model.t] = H0
	end

	j = 1
	s = zeros(Int, nspins)
	ps = zeros(Float64, nspins)
	for p in gray_code_iterator(nspins)
		s .= collect(p)
		Hj = energy(model, s)

		Pj = exp_q(-model.β * (Hj + H0), q)
		ps .= Pj * s
		Zq += Pj
		@inbounds for i in eachindex(s)
			model.x_mod[i] = model.x_mod[i] + ps[i]
		end
		i = 1
		@inbounds for k in 1:nspins-1
			for l in k+1:nspins
				model.xy_mod[i] += s[k] * ps[l]
				i += 1
			end
		end

		model.Hj_vals[j] = Hj
		model.Pj_vals[j] = Pj
		j += 1
	end

	model.Pj_vals ./= Zq
	model.x_mod ./= Zq
	model.xy_mod ./= Zq

	return nothing
end

"""
	full_measurements!(model::MaxEnt)
"""
function full_q_measurements!(model::MaxEnt)
	nspins = model.nspins
	@assert nspins < 25 "maximum full ensemble system is 24, found $nspins"

	q = model.q
	Zq = 0.0
	model.x_mod .= zeros(nspins)
	model.xy_mod .= zeros(Float64, nspins * (nspins - 1) ÷ 2)

	model.H_mean = 0.0
	H2_mean = 0.0
	model.M_mean = 0.0
	model.xyz_mod .= zeros(nspins * (nspins - 1) * (nspins - 2) ÷ 6)
	model.ones_dist_mod .= zeros(nspins + 1)

	H0 = 0.0
	if model.reg
		H0 = compute_energy_shift(model, q)
		model.H0_vals[model.t] = H0
	end
	j = 1
	s = zeros(Int, nspins)
	ps = zeros(Float64, nspins)

	for p in spin_permutations_iterator(nspins)
		s .= collect(p)
		Hj = energy(model, s)
		Pj = exp_q(-model.β * (Hj + H0), q)
		Zq += Pj
		ps .= Pj*s

		model.x_mod .= model.x_mod .+ ps
		i = 1
		for k in 1:nspins-1
			for l in k+1:nspins
				model.xy_mod[i] += s[k] * ps[l]
				i += 1
			end
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
	model.x_mod ./= Zq
	model.xy_mod ./= Zq
	model.H_mean /= Zq
	H2_mean /= Zq
	model.CV = model.β^2 * (H2_mean - model.H_mean^2)
	model.M_mean /= Zq
	model.xyz_mod ./= Zq
	model.Pj_vals ./= Zq
	model.ones_dist_mod ./= sum(model.ones_dist_mod)

	# pearson_mod!(model)
	return nothing
end
