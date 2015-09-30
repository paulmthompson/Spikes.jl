
export fit!, validate!, conmatrix

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

    Sw=Sw./k
    St=cov(m.response)
    Sb = St - Sw

    (myv, m.c.W)=eig(Sb,Sw)
    
    m.c.centroids=GroupMean*m.c.W

    nothing
    
end

function validate!{C<:LDA, V<:Training}(m::decoder{C,V})

    classes=unique(m.stimulus)
    m.predict=zeros(Float64,size(m.v.stimulus,1))

    xnew=m.v.response*m.c.W

    mydist=zeros(Float64,length(classes))
    
    for i=1:size(m.v.stimulus,1)
        for j=1:length(classes)
            mydist[j]=norm(xnew[i,:]-m.c.centroids[j,:])
        end
        m.predict[i]=classes[indmin(mydist)]
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

#=
General Methods
=#

function conmatrix{C<:classifier,V<:validation}(m::decoder{C,V},realval::Array{Float64,1})

    classes=unique(m.stimulus)
    inds=collect(1:length(classes))

    mydict=Dict{Float64,Int64}([classes[i]=>inds[i] for i=1:length(classes)])

    #y is true, x is predicted
    #sum across columns should equal 1
    conmat=zeros(Float64,length(classes),length(classes))

    for i=1:length(realval)
        yind=mydict[realval[i]]
        xind=mydict[m.predict[i]]
        conmat[yind,xind]+=1
    end

    totals_real=sum(conmat,2)
    totals_predict=sum(conmat,2)
    
    for i=1:length(classes)
        conmat[i,:]=conmat[i,:]./totals_real[i]
    end

    p_real=totals_real./sum(totals_real)
    p_predict=squeeze(totals_predict./sum(totals_predict))',2)
    
    (conmat, p_real, p_predict)
    
end
