function full_relax!(m::MaxEnt)
    mess0 = @sprintf "full_relax! "
    mess = @sprintf "%s| writing error to %s" mess0 m.err_file
    debug(LOGGER, mess)
    open(m.err_file, "w") do file
        write(file, "it,time,err_x,err_xy\n")

        err = 1.0
        t0 = time()
        told = 0

        debug(LOGGER, "aqui")
        while (m.t <= m.n_relax_steps) && (err > m.tol)
            full_iteration!(m)
            t1 = time() - t0

            # compare experimental averages (π) and model averages (Q)
            err_x = mean(abs.((m.x_obs .- m.x_mod) ./ m.x_obs))
            #err_xy = mean(abs.((m.xy_obs .- m.xy_mod) ./ m.xy_obs))
            err_xy = mean(abs.((m.pearson_obs .- m.pearson_mod) ./ m.pearson_obs))

            err = (err_x + 2 * err_xy) / 3
            write(file, "$(m.t),$(t1),$(err_x), $(err_xy)\n")
            mess1 = @sprintf "%18s|t:%5d|t(s):%6.2f" mess0 m.t t1 - told
            mess2 = @sprintf "%s|err_x:%8.5f|err_xy:%8.5f|err:%8.5f" mess1 err_x err_xy err
            info(LOGGER, mess2)

            updated_parameters!(m)
            m.t += 1
            told = t1
        end # ends relax loop

        full_iteration!(m, true)
    end # ends err file
    return nothing
end
