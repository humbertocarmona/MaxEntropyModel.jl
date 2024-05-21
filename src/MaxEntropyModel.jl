module MaxEntropyModel

using Memento
using Printf
using Graphs
using Statistics
using StatsBase: rmsd, Histogram, fit
using LinearAlgebra
using Random
using Distributions
using JSON3
using Dates
using BSON


const LOGGER = getlogger(@__MODULE__)
function __init__()
    Memento.config!("debug"; fmt="{level}: {msg}")
    # time_now = Dates.format(Dates.now(), "yy-mm-ddHH_MM_S")
    # log_file = "MaxEntropy_$(time_now).log"
    # # Create a handler for the JSON log file
    # json_handler = DefaultHandler(
    #     log_file,
    #     DictFormatter(JSON3.write)
    # )
    # # Push the json_handler to the logger
    # push!(LOGGER, json_handler)
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
    typedict
end
