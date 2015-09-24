
export pca_jpsth

function pca_jpsth(spikes_1::Array{SpikeTrain,1}, spikes_2::Array{SpikeTrain,1},trialcon::Array{Int64,1},cons::Int64,kernpca::Float64,kernpro::Float64,time::FloatRange{Float64},thispc::Int64,s::Array{Int64,1})

    binsize=time.step/time.divisor

    (mypcs, pov)=pca_session(spikes_1,kernpca,kernpca)

    println("The variance of this component is ", pov[thispc])

    myjpsth=Array(Array{Float64,3},cons)
    
    for i=1:cons
        
        myinds=find(trialcon.==i)

        thisjpsth=zeros(Float64,length(time)-1,length(time)-1,length(s))
        
        projections=project_trials(spikes_1,myinds,time,mypcs,thispc,kernpro)

        rates_spikes_2=zeros(Float64,size(projections))

        myrates=rate_KD(spikes_2,binsize,kernpro)

        count=1

        for j in s
        
            rates_spikes_2=rate_trials(myrates,myinds,time,j)

            thisjpsth[:,:,count]=cor(rates_spikes_2,projections)

            count+=1
            
        end

        myjpsth[i]=thisjpsth

    end
    
    myjpsth
    
end
