
mutable struct MaxEnt
    model::String
    runid::String

    nspins::Int64                       # number of nodes, "spins", σ≡{σ1, σ2, ... σN} 
    s::Vector{Int64}                    # store nodes states size(N), usually +1 and -1 
    Es::Float64                         # system energy corresponding to state s
    S_obs::Matrix{<:Number}                # the experimental binary matrix

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

    energy_mean::Float64
    energy_hist::Vector{Float64}
    magnetization_mean::Float64
    specific_heat::Float64

    # Random Laser specific parameters
    λwindow::Vector{<:Number}
    run_type::Char
    init_file::String
    result_file::String
    err_file::String
    comment::String
    date_today::String
    n_relax_steps::Int64
    iter_per_stage::Int64
    αh::Float64
    αJ::Float64
    γh::Float64
    γJ::Float64
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

        init_parametes!(model)
        model.β = 1.0

        model.energy_mean = 0.0 # average energy
        model.energy_hist = Float64[]
        model.magnetization_mean = 0.0 # average magnetization
        model.specific_heat = 0.0 # specific heat

        #Random Laser specific params
        model.λwindow = [1055.3, 1057.2]
        model.run_type = run_type
        model.init_file = ""
        model.comment = ""
        model.n_relax_steps = 10
        model.tol = 1.0e-4
        model.date_today = ""
        model.iter_per_stage = 1
        model.αh = 1.0
        model.αJ = 1.0
        model.γh = 0.4
        model.γJ = 0.4

        model.n_samples = 20000 * nspins
        model.n_rept = 2 * nspins
        model.n_coherence = 50 * nspins
        model.n_equilibrium = 100 * nspins

        model.result_file = "safe.json"
        model.err_file = "err.csv"
        model.mc_seed = 1234
        model.Es = energy(model)
        model.bond = make_bonds(nspins)
        model.t = 1
        return model
    end

    function MaxEnt(runid="teste", nsamples=1000, nspins=20, run_type='f')
        S = map(x -> x < 0.5 ? -1 : 1, rand(nsamples, nspins))
        nbonds = nspins * (nspins - 1) ÷ 2
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

        model.x_mod = zeros(nspins)
        model.xy_mod = zeros(nbonds)
        model.xyz_mod = zeros(size(model.xyz_obs))
        model.ones_dist_mod = zeros(size(model.ones_dist_obs))

        init_parametes!(model)
        model.β = 1.0

        model.energy_mean = 0.0 # average energy
        model.energy_hist = Float64[]
        model.magnetization_mean = 0.0 # average magnetization
        model.specific_heat = 0.0 # specific heat

        #Random Laser specific params
        model.λwindow = [1055.3, 1057.2]
        model.run_type = run_type
        model.init_file = ""
        model.comment = ""
        model.date_today = ""
        model.tol = 1.0e-4
        model.αh = 1.0
        model.αJ = 1.0
        model.γh = 0.4
        model.γJ = 0.4
        model.n_relax_steps = 20
        model.iter_per_stage = 1
        model.n_samples = 10000 * nspins
        model.n_rept = 2 * nspins
        model.n_coherence = 10 * nspins
        model.n_equilibrium = 10 * nspins
        model.result_file = "safe.json"
        model.err_file = "err.csv"
        model.bond = make_bonds(nspins)
        model.mc_seed = 1234
        model.Es = energy(model)
        model.t = 1

        return model
    end
end
