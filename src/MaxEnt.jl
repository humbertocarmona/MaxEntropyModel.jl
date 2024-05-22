mutable struct MaxEnt
	model::String
	runid::String

	nspins::Int64                       # number of nodes, "spins", σ≡{σ1, σ2, ... σN} 

	x_obs::Vector{Float64}              # observed E[ si ]
	xy_obs::Vector{Float64}             # observed E[ si⋅sj ]
	xyz_obs::Vector{Float64}            # observed E[ si⋅sj⋅sk ]
	ones_dist_obs::Vector{Float64}      # observed P[k spins up]

	x_mod::Vector{Float64}              # model computed E[si]
	xy_mod::Vector{Float64}             # model computed E[si⋅sj]
	xyz_mod::Vector{Float64}            # model computed E[si⋅sj⋅sk]
	ones_dist_mod::Vector{Float64}      # model computed P[k spins up]

	# pearson_obs::Vector{Float64}        # store observed Pearson correlation coefficient 
	# pearson_mod::Vector{Float64}        # store model Pearson correlation coefficient 

	q::Float64                          # Tsallis q
	reg::Bool                           # uses or not André's H0 regularization
	β::Float64                          # inverse temperature used to train the model
	h::Vector{Float64}                  # model local fields
	J::Vector{Float64}                  # model couplings

	Hj_vals::Vector{Float64}             # store energy values for each state
	H0_vals::Vector{Float64}            # store H0 values (Andre's correction) for each parameter set
	Pj_vals::Vector{Float64}            # store all probabilities Pj for state sj
	H_mean::Float64                     # mean system energy
	M_mean::Float64                     # system magnetization mean
	CV::Float64                         # system specific heat  β^-2 <H^2> - <H>^2

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

	"""
		MaxEnt(S::Matrix{<:Number}, runid="test", run_type='f')

		- creates a model from a matrix S
		- sets model observables
		- initialize model parameters h,J from x_obs and xy_obs

	"""
	function MaxEnt(S::Matrix{<:Number}, runid = "test", run_type = 'f')
		nspins = size(S, 2)
		model = new()
		model.model = "MaxEnt"

		model.runid = runid
		model.nspins = nspins

		model.x_obs = mean_1st_order_moments(S)
		# model.x_obs = map_to_unit_interval.(model.x_obs, -1.0, 1.0)

		model.xy_obs = mean_2nd_order_moments(S)
		# model.pearson_obs = straighten(cor(S))
		model.xyz_obs = mean_3rd_order_moments(S)
		_, model.ones_dist_obs = ones_distribution(S)

		model.x_obs[model.x_obs.==0] .+= 1e-9 # prevent inf relative error calculation
		model.xy_obs[model.xy_obs.==0] .+= 1e-9 # prevent inf relative error calculation

		model.x_mod = zeros(size(model.x_obs))
		model.xy_mod = zeros(size(model.xy_obs))
		# model.pearson_mod = zeros(size(model.pearson_obs))
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

	"""
	MaxEnt(m::MaxEnt, runid = "test", run_type = 'f')

	- just copy m and change runid and run_type

	"""
	function MaxEnt(m::MaxEnt, runid = "test", run_type = 'f')
		nspins = m.nspins
		model = copy(m)
		model.runid = runid
		model.run_type = run_type

		return model
	end

	"""
	MaxEnt(runid = "teste", nspins = 20, run_type = 'f', nsamples = 10)

	- creates a model from scratch using random S_obs

	"""
	function MaxEnt(runid = "teste", nspins = 20, run_type = 'f', nsamples = 10)
		S = map(x -> x < 0.5 ? -1 : 1, rand(nsamples, nspins))
		return MaxEnt(S, runid, run_type)
	end
end

function Base.copy(m::MaxEnt)
	model = MaxEnt(m.runid, m.nspins, m.run_type, 10)
	model.q = m.q
	model.reg = m.reg
	model.β = m.β
	model.h = copy(m.h)
	model.J = copy(m.J)
	model.Hj_vals = copy(m.Hj_vals)
	model.H0_vals = copy(m.H0_vals)
	model.Pj_vals = copy(m.Pj_vals)
	model.H_mean = m.H_mean
	model.M_mean = m.M_mean
	model.CV = m.CV
	model.run_type = m.run_type
	model.init_file = m.init_file
	model.result_file = m.result_file
	model.err_file = m.err_file
	model.comment = m.comment
	model.date_today = m.date_today
	model.λwindow = m.λwindow
	model.n_relax_steps = m.n_relax_steps
	model.ηh = m.ηh
	model.ηJ = m.ηJ
	model.γh = m.γh
	model.γJ = m.γJ
	model.α = m.α
	model.Δx = copy(m.Δx)
	model.Δxy = copy(m.Δxy)
	model.tol = m.tol
	model.n_samples = m.n_samples
	model.n_equilibrium = m.n_equilibrium
	model.n_coherence = m.n_coherence
	model.n_rept = m.n_rept
	model.mc_seed = m.mc_seed
	model.bond = copy(m.bond)
	model.t = m.t

	return model
end
