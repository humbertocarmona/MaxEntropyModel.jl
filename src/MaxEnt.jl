mutable struct MaxEnt
    model::String
    runid::String                      

    nspins::Int64                       # number of nodes, "spins", σ≡{σ1, σ2, ... σN} 
    sj::Vector{Int64}                   # store nodes states size(N), usually +1 and -1 
    Hj::Float64                         # system energy corresponding to state s
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

    q::Float64                          # Tsallis q
    reg::Bool                           # uses or not André's H0 regularization
    β::Float64                          # inverse temperature used to train the model
    h::Vector{Float64}                  # model local fields
    J::Vector{Float64}                  # model couplings

    Hj_vals::Vector{Float64}             # store energy values for each state
    H0_vals::Vector{Float64}            # store H0 values (Andre's correction) for each parameter set
    Pj_vals::Vector{Float64}            # store all probabilities Pj for state sj
    H_mean::Float64                     # mean system energy
    M_mean::Float64                         # system magnetization mean
    CV::Float64              # system specific heat  β^-2 <H^2> - <H>^2

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
        model.sj = map(x -> x < 0.5 ? 1 : -1, rand(nspins))   # random initial state
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
        model.M_mean = 0.0 # average magnetization
        model.CV = 0.0 # specific heat

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

        model.H0_vals = zeros(Float64, model.n_relax_steps)

        model.n_samples = 4000 * nspins
        model.n_rept = 1 * nspins
        model.n_coherence = 1 * nspins
        model.n_equilibrium = 100 * nspins
        model.mc_seed = 1234

        model.result_file = "safe.json"
        model.err_file = "err.csv"

        model.bond = make_bonds(nspins)
        model.Hj = energy(model)
        model.t = 1
        model.q = 1.0
        model.reg = false
        if run_type in ['f', 'q']
            model.Pj_vals = Array{Float64}(undef, 2^nspins)
            model.Hj_vals = Array{Float64}(undef, 2^nspins)
        else
            model.Pj_vals = Array{Float64}(undef, model.n_samples)
            model.Hj_vals = Array{Float64}(undef, model.n_samples)
        end


        return model
    end

    function MaxEnt(m::MaxEnt, runid="test", run_type='f')
        nspins = m.nspins
        model = MaxEnt(runid, nspins, run_type)
        model.runid = runid

        model.run_type = run_type
        model.x_obs = copy(m.x_mod)
        model.xy_obs = copy(m.xy_mod)
        model.pearson_obs = copy(m.pearson_mod)
        model.xyz_obs = copy(m.xyz_mod)
        model.ones_dist_obs = copy(m.ones_dist_mod)

        return model
    end

    function MaxEnt(runid="teste", nspins=20, run_type='f', nsamples=10)
        S = map(x -> x < 0.5 ? -1 : 1, rand(nsamples, nspins))
        return MaxEnt(S, runid, run_type)
    end
end
