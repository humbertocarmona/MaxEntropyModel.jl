module MaxEntropyModel

using Memento
using Printf
using LaTeXStrings
using Graphs
using Statistics
using StatsBase: rmsd, Histogram, fit
using LinearAlgebra
using Random
using Distributions
using JSON3
using Dates
using Makie
using BSON
using Utilities

const LOGGER = getlogger(@__MODULE__)
function __init__()
    Memento.config!("debug"; fmt="{level}: {msg}")
    # time_now = Dates.format(Dates.now(), "yy-mm-dd_HH_MM")
    # r = @sprintf "%03d" rand(collect(1:100))
    # log_file = "MaxEntropy$(r)_$(time_now).log"
    # hndlr = DefaultHandler(
    #     log_file,
    #     DefaultFormatter("{level}: {msg}")
    #     #DictFormatter(JSON3.write)
    # )
    # push!(LOGGER, hndlr)
    
    setlevel!(LOGGER, "debug")
    Memento.register(LOGGER)
end

include("MaxEnt.jl")
include("utilities.jl")
include("io_model.jl")
include("initialize_parameters.jl")
include("energy.jl")
include("full_ensemble_iteration.jl")
include("metropolis_iteration.jl")
include("update_model.jl")
include("max_entropy_relax.jl")
include("full_tsallis.jl")
include("plot_model_performance.jl")

export MaxEnt,
    init_parameters!,
    read_model,
    write_model,
    set_model!,
    full_iteration!,
    full_measurements!,
    full_q_iteration!,
    full_q_measurements!,
    metropolis_iteration!,
    metropolis_measurements!,
    max_entropy_relax!,
    straighten,
    centered_moments_obs,
    centered_moments_mod,
    energy,
    deltaEnergy, exp_q, ln_q,
    gray_code_iterator,
    plot_model_performance
end
