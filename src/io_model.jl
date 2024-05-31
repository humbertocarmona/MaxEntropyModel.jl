function save_model(model::MaxEnt, filename="none", w="bson")
    if filename == "none"
        filename = model.result_file
    end
    if w=="bson"
        filename = isnothing(findfirst(".bson", filename)) ? "$(filename).bson" : filename
        BSON.@save filename model
    else
        filename = isnothing(findfirst(".json", filename)) ? "$(filename).json" : filename
        open(filename, "w") do io
            JSON3.pretty(io, model)
        end
    end
    info(LOGGER, "wrote $filename")
end

function write_model(dic::Dict, filename::String="none";
    force=false, only=[], remove_keys=[])
    if filename == "none" && haskey(dic, :result_file)
        filename = dic[:result_file]
    end
    @assert filename != "none" "must provide a filename or key :result_file"

    filename = isnothing(findfirst(".json", filename)) ? "$(filename).json" : filename
    if size(only, 1) > 0
        push!(only, :runid)
        push!(only, :model)
        push!(only, :nspins)
        for k in keys(dic)
            if ~(k in only)
                delete!(dic, k)
            end
        end
    end
    push!(remove_keys, :Δx)
    push!(remove_keys, :Δy)
    for k in remove_keys
        delete!(dic, k)
    end

    proceed = true
    if isfile(filename) && !force
        warn(LOGGER, "$filename exists, proceed?")
        input = readline()
        proceed = (input == "y") || (input == "Y") || (input == "1")
    end
    if proceed
        open(filename, "w") do io
            JSON3.pretty(io, dic)
        end
        info(LOGGER, "wrote $filename")
    end
end

function read_model(filename::String)
    debug(LOGGER, "read_model: reading from $filename")

    obj = JSON3.read(filename)
    @assert haskey(obj, "model") "$filename must have the \"model\" key"
    @assert haskey(obj, "nspins") "$filename must have the \"nspins\" key"
    @assert obj.model == "MaxEnt"

    nspins = obj.nspins
    model = MaxEnt("teste", nspins, 'f')

    set_model!(model, Dict(obj))
    return model
end

function set_model!(model::MaxEnt, dic::Dict)
    haskey(dic, :runid) && (model.runid = dic[:runid])

    haskey(dic, :nspins) && (model.nspins = dic[:nspins])


    haskey(dic, :x_obs) && (model.x_obs = Vector{Float64}(dic[:x_obs]))
    haskey(dic, :xy_obs) && (model.xy_obs = Vector{Float64}(dic[:xy_obs]))
    haskey(dic, :xyz_obs) && (model.xyz_obs = Vector{Float64}(dic[:xyz_obs]))
    haskey(dic, :ones_dist_obs) && (model.ones_dist_obs = Vector{Float64}(dic[:ones_dist_obs]))

    haskey(dic, :x_mod) && (model.x_mod = Vector{Float64}(dic[:x_mod]))
    haskey(dic, :xy_mod) && (model.xy_mod = Vector{Float64}(dic[:xy_mod]))
    haskey(dic, :xyz_mod) && (model.xyz_mod = Vector{Float64}(dic[:xyz_mod]))
    haskey(dic, :ones_dist_mod) && (model.ones_dist_mod = Vector{Float64}(dic[:ones_dist_mod]))

    haskey(dic, :q) && (model.q = dic[:q])
    haskey(dic, :reg) && (model.reg = dic[:reg])
    haskey(dic, :β) && (model.β = dic[:β])
    haskey(dic, :h) && (model.h = Vector{Float64}(dic[:h]))
    haskey(dic, :J) && (model.J = Vector{Float64}(dic[:J]))

    haskey(dic, :Hj_vals) && (model.Hj_vals = Vector{Float64}(dic[:Hj_vals]))
    haskey(dic, :H0_vals) && (model.H0_vals = Vector{Float64}(dic[:H0_vals]))
    haskey(dic, :Pj_vals) && (model.Pj_vals = Vector{Float64}(dic[:Pj_vals]))
    haskey(dic, :H_mean) && (model.H_mean = dic[:H_mean])
    haskey(dic, :M_mean) && (model.M_mean = dic[:M_mean])
    haskey(dic, :CV) && (model.CV = dic[:CV])

    # Random Laser specific params
    if haskey(dic, :run_type)
        model.run_type = dic[:run_type][1]
    end
    haskey(dic, :init_file) && (model.init_file = dic[:init_file])
    haskey(dic, :result_file) && (model.result_file = dic[:result_file])
    haskey(dic, :err_file) && (model.err_file = dic[:err_file])
    haskey(dic, :comment) && (model.comment = dic[:comment])
    haskey(dic, :date_today) && (model.date_today = dic[:date_today])
    haskey(dic, :λwindow) && (model.λwindow = Vector{Float64}(dic[:λwindow]))
    haskey(dic, :n_relax_steps) && (model.n_relax_steps = dic[:n_relax_steps])
    haskey(dic, :ηh) && (model.ηh = dic[:ηh])
    haskey(dic, :ηJ) && (model.ηJ = dic[:ηJ])
    haskey(dic, :γh) && (model.γh = dic[:γh])
    haskey(dic, :γJ) && (model.γJ = dic[:γJ])
    haskey(dic, :α) && (model.α = dic[:α])
    # Δx::Vector{Float64}     no need
    # Δxy::Vector{Float64}    no need

    haskey(dic, :tol_x) && (model.tol_x = dic[:tol_x])
    haskey(dic, :tol_xy) && (model.tol_xy = dic[:tol_xy])
    haskey(dic, :n_samples) && (model.n_samples = dic[:n_samples])
    haskey(dic, :n_equilibrium) && (model.n_equilibrium = dic[:n_equilibrium])
    haskey(dic, :n_coherence) && (model.n_coherence = dic[:n_coherence])
    # bond no need
    haskey(dic, :t) && (model.t = dic[:t])
    haskey(dic, :n_rept) && (model.n_rept = dic[:n_rept])
    haskey(dic, :mc_seed) && (model.mc_seed = dic[:mc_seed])

    return nothing
end



