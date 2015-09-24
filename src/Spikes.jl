
module Spikes

using MAT, PyPlot, HDF5, ArrayViews, Images, DimensionalityReduction

include("types.jl")
include("getfilepath.jl")
include("behavioral.jl")
include("kinematics.jl")
include("rates.jl")
include("dimensionreduce.jl")
include("jpsth.jl")
include("neuron_dropping.jl")

end

