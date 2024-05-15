
#==
    function metropolis_iteration!(m::MaxEnt)
        Random.seed!(m.mc_seed)
        nspins = m.nspins

        m.H = energy(m)
        m.x_mod .= zeros(nspins)
        m.xy_mod .= zeros(Float64, nspins * (nspins - 1) ÷ 2)
        m.H_vals = Array{Float64}(undef, m.n_samples)

        samples = zeros(Int64, (m.n_samples, nspins))
        rept_steps = m.n_samples * m.n_coherence ÷ m.n_rept + m.n_equilibrium

        s = 1
        for _ in 1:m.n_rept
            for n in 1:rept_steps
                t = n - m.n_equilibrium
                s_flip = rand(1:m.nspins)
                flip!(m, s_flip)
                if (t > 0) && (t % m.n_coherence == 0)
                    m.x_mod .= m.x_mod .+ m.s
                    k = 1
                    for i in 1:nspins-1, j in i+1:nspins
                        m.xy_mod[k] += m.s[i] * m.s[j]
                        k += 1
                    end
                    m.H_vals[s] = m.H
                    samples[s, :] = copy(m.s)
                    s += 1
                end
            end
        end

        m.x_mod ./= m.n_samples
        m.xy_mod ./= m.n_samples

        pearson_mod!(m)

        return samples
    end
==#

#==
    function full_iteration!(m::MaxEnt)
        nspins = m.nspins
        @assert nspins < 25 "maximum full enesemble system is 24, found $nspins"

        Z = 0.0
        m.x_mod .= zeros(nspins)
        m.xy_mod .= zeros(Float64, nspins * (nspins - 1) ÷ 2)

        for p in spin_permutations_iterator(nspins)
            m.s .= collect(p)
            m.H = energy(m)
            Zx = exp(-m.β * m.H)
            Z += Zx
            m.x_mod .= m.x_mod .+ Zx .* m.s
            t = 1
            for i in 1:nspins-1, j in i+1:nspins
                m.xy_mod[t] += m.s[i] * m.s[j] * Zx
                t += 1
            end
        end
        m.x_mod ./= Z
        m.xy_mod ./= Z

        pearson_mod!(m)

        return nothing
    end
==#
