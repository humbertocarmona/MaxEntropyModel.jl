function max_entropy_relax!(model::MaxEnt)
    mess0 = @sprintf "max_entropy_relax! %s %3.2f" model.run_type model.q
    mess = @sprintf "%s| writing error to %s" mess0 model.err_file
    debug(LOGGER, mess)
    samples = nothing

    open(model.err_file, "w") do file
        write(file, "it,time,err_x,err_xy\n")

        err = 1.0
        t0 = time()
        model.t = 1
        while (model.t <= model.n_relax_steps) && (err > model.tol)
            t0 = time()
            if model.run_type == 'f'
                full_iteration!(model)
            elseif model.run_type == 'q'
                full_tsallis!(model)
            else
                metropolis_iteration!(model)
            end
            t1 = time() - t0

            # compare experimental averages (π) and model averages (Q)
            err_x = rmsd(model.x_obs, model.x_mod)
            err_xy = rmsd(model.xy_obs, model.xy_mod)

            err = (err_x + 2 * err_xy) / 3
            write(file, "$(model.t),$(t1),$(err_x),$(err_xy)\n")
            if model.t == 1 || model.t % 100 == 0
                mess1 = @sprintf "%18s|t:%5d|t(s):%6.2f" mess0 model.t t1
                mess2 = @sprintf "%s|err_x:%8.5f|err_xy:%8.5f|err:%8.5f" mess1 err_x err_xy err
                info(LOGGER, mess2)
            end
            updated_parameters!(model)
            model.t += 1
        end # ends relax loop

        if model.run_type == 'f'
            full_measurements!(model)
        elseif model.run_type == 'q'
            full_tsallis_measurements!(model)
        else
            samples = metropolis_measurements!(model)
        end
    end # ends err file

    return samples
end
