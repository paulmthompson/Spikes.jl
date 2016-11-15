
export SpikeTrain, findspikes, addevent!, rate_bin, rate_KD, Kinematic

#=
Main data container of timestamps and trial data
=#

immutable event
    inds::UnitRange{Int64}
    time::Float64
    id::UInt8
    trial::UInt16
end

type SpikeTrain
    ts::Array{Float64,1}
    trials::Array{event,1}
end

function SpikeTrain(spikes::Array{Float64,1})
    SpikeTrain(sort(spikes),Array(event,0))
end

#=
Array of spike trains coming from Array of array of time stamps
=#
function SpikeTrain(spikes::Array{Array{Float64,1},1})
    myspikes=Array(SpikeTrain,length(spikes))
    for i=1:length(spikes)
        myspikes[i]=SpikeTrain(spikes[i])
    end  
    myspikes
end

#=
Array of spike trains created from 2D array with time stamp in column 1
and neuron ID in column 2
=#
function SpikeTrain(spikes::Array{Float64,2})

    cells=unique(view(spikes,:,2))

    myspikes=[zeros(Float64,0) for i=1:length(cells)]

    myind=1

    for i=1:size(spikes,1)

        for j=1:length(cells)
            if cells[j] == spikes[i,2]
                myind=j
                break
            end
        end

        push!(myspikes[myind],spikes[i,1])

    end

    SpikeTrain(myspikes)
end

#=
Add events with events being an array of trials, with array of time stamps for each trial
=#

function addevent!(spikes::Array{SpikeTrain,1},times::Array{Array{Float64,1},1},win::Float64)

    for i=1:length(spikes)
        addevent!(spikes[i],times,win)
    end
    nothing
end

function addevent!(spikes::SpikeTrain,times::Array{Array{Float64,1},1},win::Float64)

    mysize=length(spikes.ts)
    first=searchsortedfirst(spikes.ts,times[1][1])-1
    
    if first<mysize
        for i=1:length(times)
            for j=1:length(times[i])
            
                first=searchsortedfirst(spikes.ts,(times[i][j]-win),first,mysize,Base.Forward)
                last=searchsortedfirst(spikes.ts, (times[i][j]+win),first,mysize,Base.Forward)-1

                push!(spikes.trials,event(first:last,times[i][j],j,i))       
            end
        end
    end
    nothing
end

#=
Get cell ID

(myspikes,[cell_id[round(Int64,i)] for i in cells])  

=#

#=
Rate Types
=#

abstract rate

#Bin
type rate_bin <: rate
    spikes::Array{SpikeTrain,1}
    binsize::Float64
end

#Kernel Density
type rate_KD <: rate
    spikes::Array{SpikeTrain,1}
    binsize::Float64
    kern::Float64

    function rate_KD(spikes::Array{SpikeTrain,1},binsize::Float64,kern::Float64)
        
        if kern<binsize
            kern=binsize
        end

        new(spikes,binsize,kern)
    end   
end

#=
Decoder Types
=#

export decoder, LDA, LeaveOne, Training

abstract classifier
abstract validation
abstract bias
abstract transformation

type decoder{C<:classifier,V<:validation} <: transformation
    classes::Array{Float64,1}
    c::C
    v::V
end

function decoder(c::classifier,v::validation)
    decoder(Array(Float64,0),c,v)
end

type LDA <: classifier
    myv::Array{Float64,1}
    W::Array{Float64,2}
    centroids::Array{Float64,2}
end

function LDA()
    LDA(Array(Float64,0),Array(Float64,0,0),Array(Float64,0,0))
end

type QDA <: classifier
end

type DLDA <: classifier
end

type DQDA <: classifier
end

type LeaveOne <: validation
end

type Training <: validation
    trainind::UnitRange{Int64}
    valind::UnitRange{Int64}
end

#=
Information Types
=#

export Information, QE

type Information{B<:bias,T<:transformation}
    b::B
    t::T
end

type QE <: bias
end

# Behavioral Types

type Kinematic
    px::Array{Float64,1}
    py::Array{Float64,1}
    vx::Array{Float64,1}
    vy::Array{Float64,1}
end

