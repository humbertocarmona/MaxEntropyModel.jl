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

    model.runid = obj.runid
    model.nspins = obj.nspins
    model.s = Vector{Float64}(obj.s)
    model.S_obs = reshape(S_obs, (nsamples, nspins))

    model.x_obs = Vector{Float64}(obj.x_obs)
    model.xy_obs = Vector{Float64}(obj.xy_obs)
    model.xyz_obs = Vector{Float64}(obj.xyz_obs)
    model.pearson_obs = Vector{Float64}(obj.pearson_obs)
    model.ones_dist_obs = Vector{Float64}(obj.ones_dist_obs)

    model.x_mod = Vector{Float64}(obj.x_mod)
    model.xy_mod = Vector{Float64}(obj.xy_mod)
    model.xyz_mod = Vector{Float64}(obj.xyz_mod)
    model.pearson_mod = Vector{Float64}(obj.pearson_mod)
    model.ones_dist_mod = Vector{Float64}(obj.ones_dist_mod)

    model.h = Vector{Float64}(obj.h)
    model.J = Vector{Float64}(obj.J)

    model.β = obj.β

    model.energy_mean = obj.energy_mean
    model.energy_hist = Vector{Float64}(obj.energy_hist)
    model.magnetization_mean = obj.magnetization_mean
    model.specific_heat = obj.specific_heat

    # Random Laser specific params
    model.λwindow = Vector{Float64}(obj.λwindow)
    model.run_type = obj.run_type
    model.init_file = obj.init_file
    model.result_file = obj.result_file
    model.err_file = obj.err_file
    model.comment = obj.comment
    model.date_today = obj.date_today
    model.n_relax_steps = obj.n_relax_steps
    model.tol = obj.tol
    model.iter_per_stage = obj.iter_per_stage
    model.α = obj.α
    model.γ = obj.γ
    model.n_samples = obj.n_samples
    model.n_equilibrium = obj.n_equilibrium
    model.n_coherence = obj.n_coherence
    model.bond = make_bonds(nspins)
    model.t = obj.t

    return model
end
