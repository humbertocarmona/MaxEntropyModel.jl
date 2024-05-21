struct RelaxErr
    dt::Vector{Float64}
    err_x::Vector{Float64}
    err_xy::Vector{Float64}
end


function max_entropy_relax!(model::MaxEnt)
    mess0 = @sprintf "max_entropy_relax! %s %3.2f" model.run_type model.q
    mess = @sprintf "%s| writing error to %s" mess0 model.err_file
    debug(LOGGER, mess)

    dt = Array{Float64}(undef, model.n_relax_steps)
    err_x = Array{Float64}(undef, model.n_relax_steps)
    err_xy = Array{Float64}(undef, model.n_relax_steps)
    model.H0_vals = zeros(model.n_relax_steps)


    err = 1.0
    t0 = time()
    n = 1
    while (model.t <= model.n_relax_steps) && (err > model.tol)
        t0 = time()
        if model.run_type == 'm'
            metropolis_iteration!(model)
        elseif model.q == 1.0
            full_iteration!(model)
        else
            full_q_iteration!(model)
        end

        # compare experimental averages (π) and model averages (Q)
        dt[n] = time() - t0
        err_x[n] = rmsd(model.x_obs, model.x_mod)
        err_xy[n] = rmsd(model.xy_obs, model.xy_mod)
        err = (err_x[n] + 2 * err_xy[n]) / 3

        if model.t == 1 || model.t % 10 == 0
            mess1 = @sprintf "%18s|t:%5d|dt:%6.2f" mess0 model.t dt[n]
            mess2 =
                @sprintf "%s|err_x:%8.5f|err_xy:%8.5f|err:%8.5f" mess1 err_x[n] err_xy[n] err
            info(LOGGER, mess2)
        end
        updated_parameters!(model)
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



    return RelaxErr(dt[1:n-1], err_x[1:n-1], err_xy[1:n-1])
end
