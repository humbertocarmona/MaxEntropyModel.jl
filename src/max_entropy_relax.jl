function max_entropy_relax!(model::MaxEnt; step_plot = 1)
	mess0 = @sprintf "%s|q:%6.3f" model.runid model.q
	dt = Array{Float64}(undef, model.n_relax_steps + 1)

	model.err_x = zeros(model.n_relax_steps)
	model.err_xy = zeros(model.n_relax_steps)
	model.Δx = zeros(size(model.x_obs))
	model.Δxy = zeros(size(model.xy_obs))
	model.H0_vals = zeros(model.n_relax_steps)

	t0 = time()
	n = 1
	σx_obs = std(model.x_obs) + 1.0e-12
	σxy_obs = std(model.xy_obs) + 1.0e-12
	running = true
	total_time = 0.0
	while running
		t0 = time()
		if model.run_type == 'm'
			metropolis_iteration!(model)
		elseif model.q == 1.0
			full_iteration!(model)
		else
			full_q_iteration!(model)
		end

		dt[n] = time() - t0
		total_time += dt[n]
		model.err_x[n] = rmsd(model.x_obs, model.x_mod) / σx_obs
		model.err_xy[n] = rmsd(model.xy_obs, model.xy_mod) / σxy_obs

		if model.t == 1 || model.t % step_plot == 0
			mess1 = @sprintf "%18s|t:%5d|dt:%6.3f" mess0 model.t dt[n]
			mess2 =
				@sprintf "%s|err_x:%9.6f|err_xy:%9.6f %.0f" mess1 model.err_x[n] model.err_xy[n] total_time
			info(LOGGER, mess2)

		end
		if total_time >= 1800 # save recovery every 1/2 hour
			save_model(model, "recovery_file.bson")
			total_time = 0
		end
		updated_parameters!(model)

		running = (model.err_x[n] > model.tol_x) || (model.err_xy[n] > model.tol_xy)
		running *= (model.t < model.n_relax_steps)

		model.t += 1
		n += 1
	end # ends relax loop

	if model.run_type == 'm'
		metropolis_measurements!(model)
	elseif model.q == 1.0
		full_measurements!(model)
	else
		full_q_measurements!(model)
	end
	mess3 = @sprintf "%s finished with t=%d " mess0 (model.t-1)
	info(LOGGER, mess3)

	return nothing
end
