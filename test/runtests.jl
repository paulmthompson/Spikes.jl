
module TestSetup

ENV["MPLBACKEND"] = "Agg"

using Spikes

include("rate_test.jl")

end
