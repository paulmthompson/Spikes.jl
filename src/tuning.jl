
export vec_bs, spike_count_win

#Spike Counts 

function spike_count_win(rate::SpikeTrain,t1::Float64,t2::Float64,out=zeros(Int64,length(rate.trials)))

    @inbounds for i=1:length(out)       
        mycent=rate.center[i,1]
        for j in rate.trials[i].inds
            if rate.ts[j]-mycent > t1
                if rate.ts[j]-mycent < t2
                    out[i]+=1
                end
            elseif rate.ts[j]-mycent > t2
                break
            end
        end
    end
    
    out
end

function spike_count_win(rate::SpikeTrain,t1::Array{Float64,1},t2::Array{Float64,1},out=zeros(Int64,length(rate.trials)))

    @inbounds for i=1:length(out)       
        for j in rate.trials[i].inds
            if rate.ts[j] > t1[i]
                if rate.ts[j] < t2[i]
                    out[i]+=1
                end
            elseif rate.ts[j] > t2[i]
                break
            end
        end
    end
    
    out
end

function spike_count_win_ind(rate::SpikeTrain,t1::Float64,t2::Float64,inds::Array{Int64,1},out=zeros(Int64,length(inds)))

    @inbounds for i=1:length(out)       
        mycent=rate.center[inds[i],1]
        for j in rate.trials[inds[i]].inds
            if rate.ts[j]-mycent > t1
                if rate.ts[j]-mycent < t2
                    out[i]+=1
                end
            elseif rate.ts[j]-mycent > t2
                break
            end
        end
    end
    
    out
end

function spike_count_win(rate::SpikeTrain,ts::Array{Float64,1},out=zeros(Int64,length(rate.trials),length(ts)-1))

    @inbounds for i=1:length(out)
        out[i]=0.0
    end
   @inbounds for i in 1:length(rate.trials)
        count=1
        mycent=rate.center[i,1]
        for j in rate.trials[i].inds
            myt=rate.ts[j]-mycent
            if myt>ts[count]
                k=count+1
                while k < length(ts)+1
                    if myt<ts[k]
                        out[i,k-1]+=1.0
                        count=k-1
                        break
                    end
                    k+=1
                end
            end
        end
    end
    out
end




#=
Modulation Index
Vector Method
=#

function vec_bs(rate::SpikeTrain,t1::Float64,t2::Float64,cons::Array{Int64,1},t_trials::Dict{Int64,Int64},t_dir::Array{Float64,1},inds::Array{Int64,1})
    # cons - trial type identifier for each trial of interest
    # t_trials - Maps trial type identifier to index in dictionary
    # t_dir - angle for each trial type
    # inds - maps index of t_types to index of s (first trial of interest (index 1 in t_types) may be the 4th absolute trial, so inds[1]=4)
    
    myvec1=spike_count_win(rate,t1,t2)
    
    mv1=PD_vec(myvec1,cons,t_dir,t_trials,inds)
    
    miss=0.0
    for i=1:1000
        cons2=deepcopy(cons)
        shuffle!(cons2)    
        mv2=PD_vec(myvec1,cons2,t_dir,t_trials,inds)
        if mv2>mv1
            miss+=1.0
        end
    end
    
    miss/1000
end


function PD_vec(s::Array{Int64,1},cons::Array{Int64,1},t_dir::Array{Float64,1},t_trials::Dict{Int64,Int64},inds::Array{Int64,1})
    # s - Spike counts per trial
    # cons - condition identifier for each trial of interest
    # t_dir - angle for each trial type
    # t_trials - Maps trial type identifier to index in dictionary
    # inds - maps index of t_types to index of s (first trial of interest (index 1 in t_types) may be the 4th absolute trial, so inds[1]=4)
    
    dirs=zeros(Float64,length(t_dir))
    dirn=zeros(Int64,length(t_dir))
    
    for i=1:length(cons)        
        dirs[t_trials[cons[i]]]+=s[inds[i]]    
        dirn[t_trials[cons[i]]]+=1
    end
    
    #normalize for unequal trial types
    for i=1:length(dirs)
        dirs[i]=dirs[i]/dirn[i]
    end
    
    largest_dir=maximum(dirs)
    x=0.0
    y=0.0
    
    for i=1:length(dirs)
        #dirs[i]=dirs[i]/largest_dir
        x+=cosd(t_dir[i])*dirs[i]
        y+=sind(t_dir[i])*dirs[i]
    end
     
    sqrt(x*x+y*y)
end


