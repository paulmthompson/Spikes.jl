
export SpikeTrain, findspikes, addevent!, addcenter!, rate_bin, rate_KD

#=
Main data container of timestamps and trial data
=#

immutable event
    inds::UnitRange{Int64}
    time::Float64
end

type SpikeTrain
    ts::Array{Float64,1}
    trials::Array{event,1}
    center::Array{Float64,2}
end

function SpikeTrain(spikes::Array{Float64,1},times::Array{Any,2})
    SpikeTrain(spikes,Array(event,size(times,2)),Array(Float64,size(times,2),0))
end

function findspikes(spikes::Array{Array{Float64,1},1},times::Array{Any,2},win::Float64)
    
    myspikes=Array(SpikeTrain,length(spikes))
    for i=1:length(spikes)
        myspikes[i]=SpikeTrain(spikes[i][:],times)
    end    
    
    addevent!(myspikes,times,win)

    myspikes
    
end

function findspikes(spikes::Array{Float64,2},times::Array{Any,2},win::Float64)

    spikenums=view(spikes,:,2)

    cells=unique(spikenums)
    cells=cells[2:end]
    myspikes=Array(SpikeTrain,length(cells))
    count=1
    firstind=1
    for i in cells
        lastind=searchsortedfirst(spikenums,i)
        myspikes[count]=SpikeTrain(spikes[firstind:(lastind-1),1],times)
        count+=1
        firstind=lastind
    end

    addevent!(myspikes,times,win)

    myspikes
    
end

function addevent!(spikes::Array{SpikeTrain,1},times::Array{Any,2},win::Float64)
    
    for i=1:length(spikes)
        first=searchsortedfirst(spikes[i].ts,times[1][1])-1
        mysize=length(spikes[i].ts)
        if first<mysize
            for j=1:size(times,2)
                first=searchsortedfirst(spikes[i].ts,(times[j][1]-win),first,mysize,Base.Forward)
                last=searchsortedfirst(spikes[i].ts, (times[j][end]+win),first,mysize,Base.Forward)-1
                spikes[i].trials[j]=Spikes.event(first:last,times[j][1])
            end
        end
    end

    nothing
end

function addcenter!(spikes::Array{SpikeTrain,1},center::Array{Float64,1})
      
    for i=1:length(spikes)
        spikes[i].center=hcat(spikes[i].center,center)
    end

    nothing
end

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

type LDA <: classifier
    W::Array{Float64,2}
    centroids::Array{Float64,2}
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

type Information{B<:bias,T<:transformation}
    b::B
    t::T
end

type QE <: bias
end


    

