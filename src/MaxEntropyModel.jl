module MaxEntropyModel

using Memento
using Printf
using Graphs
using Statistics
using StatsBase: Histogram, fit
using LinearAlgebra
using Random
using Distributions
using JSON3
using Dates

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
include("energy.jl")
include("full_ensemble_iteration.jl")
include("update_model.jl")
include("full_ensemble_relax.jl")

export MaxEnt,
    read_model,
    write_model,
    full_relax!,
    full_iteration!,
    straighten,
    centered_moments_obs,
    centered_moments_mod

end
