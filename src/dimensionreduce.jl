
export pca_session

function pca_session(spikes::Array{SpikeTrain,1},binsize::Float64,kern::Float64)

    myrate=rate_KD(spikes,binsize,kern)
    
    pcamat=zeros(Float64,length(spikes[1].trials[1].time:binsize:spikes[1].trials[end].time)-1,length(spikes))

for i=1:length(spikes)
    pcamat[:,i]=rate_session(myrate,i)
end

    mypca=pca(pcamat)

    (mypca.rotation, mypca.proportion_of_variance)
    
end

function project_trials(spikes::Array{SpikeTrain,1},inds::Array{Int64,1},timeseq::FloatRange{Float64},components::Array{Float64,2},mypc::Int64,kern::Float64)

    binsize=timeseq.step/timeseq.divisor
    
    myrate=rate_KD(spikes,binsize,kern)
    
    pca_vec=components[:,mypc]
    
    projections=zeros(Float64,length(inds),length(timeseq)-1)

    trialrates=zeros(Float64,length(timeseq)-1,length(spikes))
    
    for i=1:length(inds)

        trialrates=rate_window_pop(myrate,inds[i],timeseq)
        
        projections[i,:]=trialrates*pca_vec
    
    end

    projections
    
end
