function max_entropy_relax!(m::MaxEnt)
    mess0 = @sprintf "max_entropy_relax! %s " m.run_type
    mess = @sprintf "%s| writing error to %s" mess0 m.err_file
    debug(LOGGER, mess)
    samples = nothing

    open(m.err_file, "w") do file
        write(file, "it,time,err_x,err_xy\n")

        err = 1.0
        t0 = time()
        m.t = 1
        while (m.t <= m.n_relax_steps) && (err > m.tol)
            t0 = time()
            if m.run_type == 'f'
                full_iteration!(m)
            elseif m.run_type == 'q'
                full_tsallis!(m)
            else
                metropolis_iteration!(m)
            end
            t1 = time() - t0

            # compare experimental averages (π) and model averages (Q)
            err_x = rmsd(m.x_obs, m.x_mod)
            err_xy = rmsd(m.xy_obs, m.xy_mod)

            err = (err_x + 2 * err_xy) / 3
            write(file, "$(m.t),$(t1),$(err_x),$(err_xy)\n")
            if m.t == 1 || m.t % 100 == 0
                mess1 = @sprintf "%18s|t:%5d|t(s):%6.2f" mess0 m.t t1
                mess2 = @sprintf "%s|err_x:%8.5f|err_xy:%8.5f|err:%8.5f" mess1 err_x err_xy err
                info(LOGGER, mess2)
            end
            updated_parameters!(m)
            m.t += 1
        end # ends relax loop

        if m.run_type == 'f'
            full_measurements!(m)
        elseif m.run_type == 'q'
            full_tsallis_measurements!(m)
        else
            samples = metropolis_measuremments!(m)
        end
    end # ends err file

    return samples
end
