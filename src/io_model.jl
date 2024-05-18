using JSON3
function write_model(model::MaxEnt, filename::String = "none"; force = false)
	if filename == "none"
		filename = model.result_file
	end
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

function write_model(dic::Dict, filename::String = "none";
	force = false, only = [], remove_keys=[])
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
	push!(remove_keys,:Δx)
	push!(remove_keys,:Δy)
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
	model = MaxEnt("teste", nspins,'f')

	haskey(obj, "runid") && (model.runid = obj.runid)
	haskey(obj, "nspins") && (model.nspins = obj.nspins)
	haskey(obj, "s") && (model.sj = Vector{Float64}(obj.sj))

	if haskey(obj, "S_obs")
		S_obs = Vector{Int64}(obj.S_obs)
		nsamples = length(S_obs) ÷ nspins
		model.S_obs = reshape(S_obs, (nsamples, nspins))
	end

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

function set_model!(model::MaxEnt, dic::Dict)
	haskey(dic, :runid) && (model.runid = dic[:runid])

	haskey(dic, :nspins) && (model.nspins = dic[:nspins])
	haskey(dic, :sj) && (model.sj = Vector{Float64}(dic[:sj]))
	haskey(dic, :Hj) && (model.Hj = Vector{Float64}(dic[:Hj]))
	if haskey(dic, :S_obs)
		S_obs = Vector{Int64}(obj[:S_obs])
		nsamples = length(S_obs) ÷ nspins
		model.S_obs = reshape(S_obs, (nsamples, nspins))
	end

	haskey(dic, :x_obs) && (model.x_obs = Vector{Float64}(dic[:x_obs]))
	haskey(dic, :xy_obs) && (model.xy_obs = Vector{Float64}(dic[:xy_obs]))
	haskey(dic, :xyz_obs) && (model.xyz_obs = Vector{Float64}(dic[:xyz_obs]))
	haskey(dic, :pearson_obs) && (model.pearson_obs = Vector{Float64}(dic[:pearson_obs]))
	haskey(dic, :ones_dist_obs) && (model.ones_dist_obs = Vector{Float64}(dic[:ones_dist_obs]))

	haskey(dic, :x_mod) && (model.x_mod = Vector{Float64}(dic[:x_mod]))
	haskey(dic, :xy_mod) && (model.xy_mod = Vector{Float64}(dic[:xy_mod]))
	haskey(dic, :xyz_mod) && (model.xyz_mod = Vector{Float64}(dic[:xyz_mod]))
	haskey(dic, :pearson_mod) && (model.pearson_mod = Vector{Float64}(dic[:pearson_mod]))
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

	haskey(dic, :tol) && (model.tol = dic[:tol])
	haskey(dic, :n_samples) && (model.n_samples = dic[:n_samples])
	haskey(dic, :n_equilibrium) && (model.n_equilibrium = dic[:n_equilibrium])
	haskey(dic, :n_coherence) && (model.n_coherence = dic[:n_coherence])
	# bond no need
	haskey(dic, :t) && (model.t = dic[:t])
	haskey(dic, :n_rept) && (model.n_rept = dic[:n_rept])
	haskey(dic, :mc_seed) && (model.mc_seed = dic[:mc_seed])

	return nothing
end

function set_model!(model::MaxEnt, json_string::String)
	obj = JSON3.read(json_string)
	set_model!(model, Dict(obj))

	return nothing
end

function set_model!(input_file::String)
	obj = JSON3.read(filename)
	@assert haskey(obj, "model") "$filename must have the \"model\" key"
	@assert haskey(obj, "nspins") "$filename must have the \"nspins\" key"
	@assert obj.model == "MaxEnt"
	model = MaxEnt("teste", obj["nspins"],'f')
	set_model!(model, Dict(obj))

	return model
end

@inline function typedict(x::T) where {T} 
	return  Dict(fn=>getfield(x, fn) for fn ∈ fieldnames(T))
end

