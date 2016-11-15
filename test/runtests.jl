
module TestSetup

ENV["MPLBACKEND"] = "Agg"

using Spikes

include("type_test.jl")
include("rate_test.jl")

end
