
module Spikes

using MAT, PyPlot, HDF5, ArrayViews, Images, DimensionalityReduction, StatsBase, PyCall

include("types.jl")
include("rates.jl")
include("dimensionreduce.jl")
include("jpsth.jl")
include("neuron_dropping.jl")
include("information.jl")
include("classifiers.jl")

end

