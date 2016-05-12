
module Spikes

using MAT, PyPlot, HDF5, ArrayViews, Images, DimensionalityReduction, StatsBase, PyCall, HypothesisTests

@pyimport matplotlib.collections as C

include("types.jl")
include("rates.jl")
include("tuning.jl")
include("dimensionreduce.jl")
include("jpsth.jl")
include("neuron_dropping.jl")
include("information.jl")
include("classifiers.jl")
include("kinematics.jl")

end

