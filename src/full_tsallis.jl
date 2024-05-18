
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

	# ---- using deltaEnergy
	# sj = ones(Int, nspins)
	# model.sj .= copy(sj)
	# model.Hj = energy(model)

	j = 1
	for p in gray_code_iterator(nspins)
		model.sj .= collect(p)
		model.Hj = energy(model)
		
		# ---- using deltaEnergy
		# sj .= collect(p)
		# i = findfirst(x->x!=0, sj - model.sj)
		# if ~isnothing(i)
		# 	model.Hj += deltaEnergy(model,i)
		# 	model.sj .= sj
		# end

		Pj = exp_q(-model.β * (model.Hj + H0), q)
		Zq += Pj
		model.x_mod .= model.x_mod .+ Pj .* model.sj
		i = 1
		@simd for k in 1:nspins-1
			@simd for l in k+1:nspins
				@inbounds model.xy_mod[i] += model.sj[k] * model.sj[l] * Pj
				i += 1
			end
		end

		model.Hj_vals[j] = model.Hj
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

TBW
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

	H0 = compute_energy_shift(model, q)
	model.H0_vals[1] = H0
	j = 1
	for s in spin_permutations_iterator(nspins)
		model.sj .= collect(s)
		model.Hj = energy(model)
		if model.reg
			Pj = exp_q(-model.β * (model.Hj + H0), q)
		else
			Pj = exp_q(-model.β * model.Hj, q)
		end
		Zq += Pj

		model.x_mod .= model.x_mod .+ Pj .* model.sj
		i = 1
		for k in 1:nspins-1, l in k+1:nspins
			model.xy_mod[i] += model.sj[k] * model.sj[l] * Pj
			i += 1
		end

		model.H_mean += model.Hj * Pj
		H2_mean += model.Hj * model.Hj * Pj
		model.M_mean += sum(model.sj) * Pj
		i = 1
		for k in 1:nspins-2
			for l in k+1:nspins-1
				for m in l+1:nspins
					model.xyz_mod[i] += model.sj[k] * model.sj[l] * model.sj[m] * Pj
					i += 1
				end
			end
		end
		k = count(isone.(model.sj))
		model.ones_dist_mod[k+1] += Pj

		model.Hj_vals[j] = model.Hj
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

	pearson_mod!(model)
	return nothing
end
