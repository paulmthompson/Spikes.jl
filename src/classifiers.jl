
export fit!, validate!

#=
Maximum Liklihood

Linear Discriminant Analysis
Quadratic Discriminant Analysis
Diagonal Discriminant Analysis
Diagonal Quadratic Discriminant Analysis
Shrinkage?
=#

#=
LDA
=#

function fit!{C<:LDA, V<:validation}(m::decoder{C,V})

    classes=unique(m.stimulus)
    k=length(classes)

    nGroup=zeros(Int64,k)
    GroupMean=zeros(Float64,k,size(m.response,2))
    Sw=zeros(Float64,size(m.response,2),size(m.response,2))
    m.c.W=zeros(Float64,k,size(m.response,2)+1)

    for i=1:k
        Group = m.stimulus.==classes[i]
        nGroup[i]=sum(Group)

        GroupMean[i,:]=mean(m.response[Group,:],1)

        Sw += cov(m.response[Group,:])
        
    end

    Sw=Sw.e/k
    St=cov(m.response)
    Sb = St - Sw

    (myv, m.c.W)=eig(Sb,Sw)
    
    PriorProb = nGroup / size(m.response,1)

    for i = 1:k

        Temp = GroupMean[i,:] * m.c.W * m.c.W'
    
        m.c.W[i,1] = (-0.5 * Temp * GroupMean[i,:]')[1] + log(PriorProb[i])
    
        m.c.W[i,2:end] = Temp
        
    end

    nothing
    
end

function validate!{C<:LDA, V<:Training}(m::decoder{C,V})

    classes=unique(m.stimulus)
    m.predict=zeros(Float64,size(m.v.stimulus,1))
    
    for i=1:size(m.v.stimulus,1)
        project=m.c.w*m.v.response[i,:]
        m.predict=classes[minind(abs(project-classes))]
    end

    nothing
    
end

#=
Nearest Neighbor Classification

Quian Quiroga et al 2006
=#

#=
Naive Bayesian
=#
