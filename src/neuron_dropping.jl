
export neuron_dropping, neuron_linear_regression, population_linear_regression

function neuron_dropping(spikes::Array{SpikeTrain,1},analog::Array{Float64,2},binsize::Float64,nlags::UnitRange{Int64},firstwin::UnitRange{Int64},secondwin::UnitRange{Int64},orderinds::Array{Int64,1})

    neuroncorr=neuron_linear_regression(spikes,analog,binsize,nlags,firstwin)

    # Y is first x minutes of analog data 
    Y = analog[firstwin, :]
    mywin1=collect(linspace((firstwin[1]-nlags[end]-1)*binsize, (firstwin[end]+1)*binsize, length(firstwin)+length(nlags)))

    # Y2 is second x minutes of analog data 
    Y2 = analog[secondwin, :]
    mywin2=collect(linspace((secondwin[1]-nlags[end]-1)*binsize, (secondwin[end]+1)*binsize, length(secondwin)+length(nlags)))
    
    #find order of neurons from highest correlation to lowest correlation
    myinds=sortperm(squeeze(maximum(neuroncorr[:,orderinds],2),2), rev=true)
    
    #setup master neuron firing rate matrix, which is ordered X1, X2, X3 like above (so columns = neuron * (nlags+1), with order defined by highest correlation
    #first column and last column are ones

    neurbin1=Array(Float64,size(Y,1),(length(nlags)*length(spikes))+2)
    neurbin2=Array(Float64,size(Y,1),(length(nlags)*length(spikes))+2)

    neurbin1[:,1]=ones(Float64,size(neurbin1,1))
    neurbin2[:,1]=ones(Float64,size(neurbin2,1))
    
    count=2
    
    for i in myinds

        ts1=hist(spikes[i].ts,mywin1)[2]
        ts2=hist(spikes[i].ts,mywin2)[2]
        
        for j in nlags
            neurbin1[:,count]=ts1[(nlags[end]-j+1):(end-j)]
            neurbin2[:,count]=ts2[(nlags[end]-j+1):(end-j)]
            count+=1
        end

    end

    neurbin1[:,end]=ones(Float64,size(neurbin1,1))
    neurbin2[:,end]=ones(Float64,size(neurbin2,1))

    mincorr=zeros(Float64,length(spikes),size(analog,2))

    for i=length(spikes):-1:1
        myrange=((i-1)*(length(nlags))+2):size(neurbin1,2)
        X=view(neurbin1,:,myrange)
        coeffs= X \ Y
        X2=view(neurbin2,:,myrange)
        Z= X2 * coeffs
        mincorr[length(spikes)-i+1,:]=abs(diag(cor(Z,Y2)))
    end
      
    #max corr iteration

    maxcorr=zeros(Float64,length(spikes),size(analog,2))

    for i=1:length(spikes)

        myrange=1:((i*length(nlags))+1)
        X=view(neurbin1,:,myrange)
        coeffs= X \ Y
        X2=view(neurbin2,:,myrange)
        Z = X2 * coeffs
        maxcorr[i,:]=abs(diag(cor(Z,Y2)))
        
    end

    hcat(mincorr,maxcorr)
end

function neuron_linear_regression(spikes::Array{SpikeTrain,1},analog::Array{Float64,2},binsize::Float64,nlags::UnitRange{Int64},firstwin::UnitRange{Int64})

    neuroncorr=Array(Float64,length(spikes),size(analog,2))
  
    # Y is first x minutes of analog data 
    Y = analog[firstwin, :]
    mywin1=collect(linspace((firstwin[1]-nlags[end]-1)*binsize, (firstwin[end]+1)*binsize, length(firstwin)+length(nlags)))

    X=zeros(Float64,size(Y,1),length(nlags)+1)
    X[:,1]=ones(Float64,size(Y,1),1)
    
    for i=1:length(spikes)
    
        #bin time stamps for first part
        ts=hist(spikes[i].ts,mywin1)[2]
        count=2
        for j in nlags
            X[:,count]=ts[(nlags[end]-j+1):(end-j)]
            count+=1
        end

        coeffs = X \ Y

        Z = X * coeffs

        neuroncorr[i,:]=abs(diag(cor(Z,Y)))

    end

    neuroncorr

end

function population_linear_regression(spikes::Array{SpikeTrain,1},analog::Array{Float64,2},binsize::Float64,nlags::UnitRange{Int64},firstwin::UnitRange{Int64})
  
    # Y is first x minutes of analog data 
    Y = analog[firstwin, :]
    mywin1=collect(linspace((firstwin[1]-nlags[end]-1)*binsize, (firstwin[end]+1)*binsize, length(firstwin)+length(nlags)))

    X=zeros(Float64,size(Y,1),(length(nlags)*length(spikes))+1)
    X[:,1]=ones(Float64,size(Y,1),1)

    count=2
    
    for i=1:length(spikes)

        ts1=hist(spikes[i].ts,mywin1)[2]
        
        for j in nlags
            X[:,count]=ts1[(nlags[end]-j+1):(end-j)]
            count+=1
        end

    end

    coeffs = X \ Y

    Z = X * coeffs

    neuroncorr=abs(diag(cor(Z,Y)))

end

function population_linear_regression_predict(spikes::Array{SpikeTrain,1},analog::Array{Float64,2},binsize::Float64,nlags::UnitRange{Int64},firstwin::UnitRange{Int64},secondwin::UnitRange{Int64})

        # Y is first x minutes of analog data 
    Y = analog[firstwin, :]
    mywin1=collect(linspace((firstwin[1]-nlags[end]-1)*binsize, (firstwin[end]+1)*binsize, length(firstwin)+length(nlags)))
    Y2 = analog[secondwin,:]

    X=zeros(Float64,size(Y,1),(length(nlags)*length(spikes))+1)
    X[:,1]=ones(Float64,size(Y,1),1)

    count=2
    
    for i=1:length(spikes)

        ts1=hist(spikes[i].ts,mywin1)[2]
        
        for j in nlags
            X[:,count]=ts1[(nlags[end]-j+1):(end-j)]
            count+=1
        end

    end

    coeffs = X \ Y

    X2 = zeros(Float64,size(Y2,1),(length(nlags)*length(spikes))+1)
    X2[:,1]=ones(Float64,size(Y2,1),1)

    mywin2=collect(linspace((secondwin[1]-nlags[end]-1)*binsize, (secondwin[end]+1)*binsize, length(secondwin)+length(nlags)))

    count=2
    
    for i=1:length(spikes)

        ts2=hist(spikes[i].ts,mywin2)[2]
        
        for j in nlags
            X2[:,count]=ts2[(nlags[end]-j+1):(end-j)]
            count+=1
        end

    end

    Z = X2 * coeffs

    neuroncorr=abs(diag(cor(Z,Y2)))

end
