
#=
Maximum Liklihood

Linear Discriminant Analysis
Quadratic Discriminant Analysis
Diagonal Discriminant Analysis
Diagonal Quadratic Discriminant Analysis
=#

#=
LDA
=#

function LDA(response::Array{Float64,2},stimulus::Array{Float64,1})

    classes=unique(stimulus)
    k=length(classes)

    nGroup=zeros(Int64,k)
    GroupMean=zeros(Float64,k,size(response,2))
    PooledCov=zeros(Float64,size(response,2),size(response,2))
    W=zeros(k,size(response,2)+1)

    for i=1:k
        Group = stimulus==classes[i]
        nGroup[i]=sum(Group)

        GroupMean[i,:]=mean(response[Group,:])

        PooledCov= PooledCov + ((nGroup[i] - 1) / (size(response,1) - k) ) .* cov(response[Group,:])

    end

    PriorProb = nGroup / size(response,1)

    for i = 1:k,

        Temp = GroupMeanp[i,:] / PooledCov
    
        W[i,1] = -0.5 * Temp * GroupMean[i,:]' + log(PriorProb[i])
    
        W[i,2:end] = Temp
        
    end

    W
    
end


#=
Nearest Neighbor Classification

Quian Quiroga et al 2006
=#

#=
Naive Bayesian
=#
