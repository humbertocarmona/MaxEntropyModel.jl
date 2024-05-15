
mutable struct MaxEnt
    model::String
    runid::String

    nspins::Int64                       # number of nodes, "spins", σ≡{σ1, σ2, ... σN} 
    s::Vector{Int64}                    # store nodes states size(N), usually +1 and -1 
    H::Float64                          # system energy corresponding to state s
    S_obs::Matrix{<:Number}             # the experimental binary matrix

    x_obs::Vector{Float64}              # observed E[ si ]
    xy_obs::Vector{Float64}             # observed E[ si⋅sj ]
    xyz_obs::Vector{Float64}            # observed E[ si⋅sj⋅sk ]
    pearson_obs::Vector{Float64}        # store observed Pearson correlation coefficient 
    ones_dist_obs::Vector{Float64}      # observed P[k spins up]

    x_mod::Vector{Float64}              # model computed E[si]
    xy_mod::Vector{Float64}             # model computed E[si⋅sj]
    xyz_mod::Vector{Float64}            # model computed E[si⋅sj⋅sk]
    pearson_mod::Vector{Float64}        # store model Pearson correlation coefficient 
    ones_dist_mod::Vector{Float64}      # model computed P[k spins up]

    h::Vector{Float64}                  # model local fields
    J::Vector{Float64}                  # model couplings

    β::Float64                          # inverse temperature used to train the model

    H_mean::Float64                     # mean system energy
    H_vals::Vector{Float64}             # store energy values for each state
    H0_vals::Vector{Float64}            # store H0 values (Andre's correction) for each parameter set

    PE_weights::Vector{Float64}                 # store density of states P(E)
    PE_edges::Vector{Float64}
    Hmin::Float64                       # minimum energy for the density of states
    Hmax::Float64                       # maximum energy for the density of states

    magnetization_mean::Float64         # system magnetization mean
    specific_heat::Float64              # system specific heat  β^-2 <H^2> - <H>^2

    # Random Laser specific parameters
    run_type::Char
    init_file::String
    result_file::String
    err_file::String
    comment::String
    date_today::String
    λwindow::Vector{<:Number}
    n_relax_steps::Int64
    ηh::Float64
    ηJ::Float64
    γh::Float64
    γJ::Float64
    α::Float64              # update inertia
    Δx::Vector{Float64}     # ..
    Δxy::Vector{Float64}    # ..

    tol::Float64
    n_samples::Int64
    n_equilibrium::Int64
    n_coherence::Int64
    n_rept::Int64
    mc_seed::Int64
    # helper bond[i,j] = t, the t-th element in J or xy
    bond::Matrix{Int64}
    t::Int64

    function MaxEnt(S::Matrix{<:Number}, runid="test", run_type='f')
        nspins = size(S, 2)
        model = new()
        model.model = "MaxEnt"

        model.runid = runid
        model.nspins = nspins
        model.s = map(x -> x < 0.5 ? 1 : -1, rand(nspins))   # random initial state
        model.S_obs = copy(S)

        model.x_obs = mean_1st_order_moments(S)
        # model.x_obs = map_to_unit_interval.(model.x_obs, -1.0, 1.0)

        model.xy_obs = mean_2nd_order_moments(S)
        model.pearson_obs = straighten(cor(S))
        model.xyz_obs = mean_3rd_order_moments(S)
        _, model.ones_dist_obs = ones_distribution(S)

        model.x_obs[model.x_obs.==0] .+= 1e-9 # prevent inf relative error calculation
        model.xy_obs[model.xy_obs.==0] .+= 1e-9 # prevent inf relative error calculation


        model.x_mod = zeros(size(model.x_obs))
        model.xy_mod = zeros(size(model.xy_obs))
        model.pearson_mod = zeros(size(model.pearson_obs))
        model.xyz_mod = zeros(size(model.xyz_obs))
        model.ones_dist_mod = zeros(size(model.ones_dist_obs))

        model.Δx = zeros(size(model.x_obs))
        model.Δxy = zeros(size(model.xy_obs))
        init_parameters!(model)
        model.β = 1.0

        model.H_mean = 0.0 # average energy
        model.H_vals = Float64[]
        model.magnetization_mean = 0.0 # average magnetization
        model.specific_heat = 0.0 # specific heat

        #Random Laser specific params
        model.λwindow = [1055.3, 1057.2]
        model.run_type = run_type
        model.init_file = ""
        model.comment = ""
        model.date_today = ""

        model.ηh = 1.0
        model.ηJ = 1.1
        model.γh = 0.2
        model.γJ = 0.2
        model.α = 0.1
        model.Δx = zeros(size(model.x_obs))
        model.Δxy = zeros(size(model.xy_obs))
        model.n_relax_steps = 50
        model.tol = 1.0e-6

        model.H0_vals = Array{Float64}(undef, model.n_relax_steps)

        model.n_samples = 4000 * nspins
        model.n_rept = 1 * nspins
        model.n_coherence = 1 * nspins
        debug(LOGGER, "n_coherene=$(model.n_coherence)")
        model.n_equilibrium = 100 * nspins
        model.mc_seed = 1234

        model.result_file = "safe.json"
        model.err_file = "err.csv"

        model.bond = make_bonds(nspins)
        model.H = energy(model)
        model.t = 1

        model.Hmin = -20.0
        model.Hmax = 20.0
        nbins = 2 * round(Int64, 2.0^(model.nspins / 3)) # Rice rule 
        model.PE_edges = LinRange(model.Hmin, model.Hmax, nbins + 1)
        model.PE_weights = zeros(Float64, nbins)
        return model
    end

    function MaxEnt(runid="teste", nsamples=1000, nspins=20, run_type='f')
        S = map(x -> x < 0.5 ? -1 : 1, rand(nsamples, nspins))
        return MaxEnt(S, runid, run_type)
    end
end
