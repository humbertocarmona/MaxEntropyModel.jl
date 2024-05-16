using JSON3
function write_model(model::MaxEnt, filename::String; force=false)
    filename = isnothing(findfirst(".json", filename)) ? "$(filename).json" : filename

    proceed = true
    if isfile(filename) && !force
        warn(LOGGER, "$filename exists, proceed?")
        input = readline()
        proceed = (input == "y") || (input == "Y") || (input == "1")
    end
    if proceed
        open(filename, "w") do io
            JSON3.pretty(io, model)
        end
        info(LOGGER, "wrote $filename")
    end
end

function read_model(filename::String)
    debug(LOGGER, "read_model: reading from $filename")

    obj = JSON3.read(filename)
    @assert haskey(obj, "model") "$filename must have the \"model\" key"
    @assert haskey(obj, "nspins") "$filename must have the \"nspins\" key"
    @assert haskey(obj, "S_obs") "$filename must have the \"S_obs\" key"
    @assert obj.model == "MaxEnt"

    nspins = obj.nspins
    model = MaxEnt()
    S_obs = Vector{Int64}(obj.S_obs)
    nsamples = length(S_obs) ÷ nspins

    haskey(obj, "runid") && (model.runid = obj.runid)
    haskey(obj, "nspins") && (model.nspins = obj.nspins)
    haskey(obj, "s") && (model.sj = Vector{Float64}(obj.sj))
    haskey(obj, "S_obs") && (model.S_obs = reshape(S_obs, (nsamples, nspins)))

    haskey(obj, "x_obs") && (model.x_obs = Vector{Float64}(obj.x_obs))
    haskey(obj, "xy_obs") && (model.xy_obs = Vector{Float64}(obj.xy_obs))
    haskey(obj, "xyz_obs") && (model.xyz_obs = Vector{Float64}(obj.xyz_obs))
    haskey(obj, "pearson_obs") && (model.pearson_obs = Vector{Float64}(obj.pearson_obs))
    haskey(obj, "ones_dist_obs") && (model.ones_dist_obs = Vector{Float64}(obj.ones_dist_obs))

    haskey(obj, "x_mod") && (model.x_mod = Vector{Float64}(obj.x_mod))
    haskey(obj, "xy_mod") && (model.xy_mod = Vector{Float64}(obj.xy_mod))
    haskey(obj, "xyz_mod") && (model.xyz_mod = Vector{Float64}(obj.xyz_mod))
    haskey(obj, "pearson_mod") && (model.pearson_mod = Vector{Float64}(obj.pearson_mod))
    haskey(obj, "ones_dist_mod") && (model.ones_dist_mod = Vector{Float64}(obj.ones_dist_mod))

    haskey(obj, "h") && (model.h = Vector{Float64}(obj.h))
    haskey(obj, "J") && (model.J = Vector{Float64}(obj.J))

    haskey(obj, "β") && (model.β = obj.β)
    haskey(obj, "q") && (model.q = obj.q)
    haskey(obj, "reg") && (model.reg = obj.reg)

    haskey(obj, "H_mean") && (model.H_mean = obj.H_mean)
    haskey(obj, "Hj_vals") && (model.Hj_vals = Vector{Float64}(obj.Hj_vals))
    haskey(obj, "Pj_vals") && (model.Pj_vals = Vector{Float64}(obj.Pj_vals))
    haskey(obj, "H0_vals") && (model.H0_vals = Vector{Float64}(obj.H0_vals))
    haskey(obj, "M_mean") && (model.M_mean = obj.M_mean)
    haskey(obj, "CV") && (model.CV = obj.CV)

    # Random Laser specific params
    haskey(obj, "λwindow") && (model.λwindow = Vector{Float64}(obj.λwindow))
    if haskey(obj, "run_type")
        model.run_type = obj.run_type[1]
    end
    haskey(obj, "init_file") && (model.init_file = obj.init_file)
    haskey(obj, "result_file") && (model.result_file = obj.result_file)
    haskey(obj, "err_file") && (model.err_file = obj.err_file)
    haskey(obj, "comment") && (model.comment = obj.comment)
    haskey(obj, "date_today") && (model.date_today = obj.date_today)
    haskey(obj, "n_relax_steps") && (model.n_relax_steps = obj.n_relax_steps)
    haskey(obj, "tol") && (model.tol = obj.tol)
    haskey(obj, "ηh") && (model.ηh = obj.ηh)
    haskey(obj, "γh") && (model.γh = obj.γh)
    haskey(obj, "ηJ") && (model.ηJ = obj.ηh)
    haskey(obj, "γJ") && (model.γJ = obj.γh)
    haskey(obj, "ηJ") && (model.ηJ = obj.ηJ)
    haskey(obj, "γJ") && (model.γJ = obj.γJ)
    haskey(obj, "α") && (model.α = obj.α)
    haskey(obj, "n_samples") && (model.n_samples = obj.n_samples)
    haskey(obj, "n_equilibrium") && (model.n_equilibrium = obj.n_equilibrium)
    haskey(obj, "n_coherence") && (model.n_coherence = obj.n_coherence)
    haskey(obj, "bond") && (model.bond = make_bonds(nspins))
    haskey(obj, "t") && (model.t = obj.t)
    haskey(obj, "mc_seed") && (model.mc_seed = obj.mc_seed)
    haskey(obj, "n_rept") && (model.n_rept = obj.n_rept)

    model.Hj = energy(model)
    return model
end
