
export fit!, project, conmatrix

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

function fit!{C<:LDA, V<:validation}(m::decoder{C,V},response::Array{Float64,2},stimulus::Array{Float64,1})

    m.classes=unique(stimulus)
    k=length(m.classes)

    #if any stimulus appears <10 times, ignore it
    keep=Array(Int64,0)
    for i=1:k
        Group = sum(stimulus.==m.classes[i])
        if Group<10
        else
            push!(keep,i)
        end
    end

    m.classes=m.classes[keep]
    k=length(m.classes)
        

    nGroup=zeros(Int64,k)
    GroupMean=zeros(Float64,k,size(response,2))
    Sw=zeros(Float64,size(response,2),size(response,2))
    m.c.W=zeros(Float64,k,size(response,2))
    covmat=zeros(Float64,size(response,2),size(response,2),k)

    for i=1:k
        Group = stimulus.==m.classes[i]
        nGroup[i]=sum(Group)

        GroupMean[i,:]=mean(response[Group,:],1)

        covmat[:,:,i] = cov(response[Group,:])
        
    end

    for i=1:k
        Sw+=(nGroup[i]/sum(nGroup)).*covmat[:,:,i]
    end
    
    St=cov(response)
    Sb = St - Sw

    (m.c.myv, W)=eig(Sb,Sw)

    mostvar=find((m.c.myv./sum(m.c.myv)).>.01)

    m.c.W=W[:,mostvar]
    
    m.c.centroids=GroupMean*m.c.W

    nothing
    
end

function project{C<:LDA, V<:validation}(m::decoder{C,V},response::Array{Float64,2},stimulus::Array{Float64,1})

    predict=zeros(Float64,size(stimulus,1))

    xnew=response*m.c.W

    mydist=zeros(Float64,length(m.classes))
    
    for i=1:size(stimulus,1)
        for j=1:length(m.classes)
            mydist[j]=norm(xnew[i,:]-m.c.centroids[j,:])
        end
        predict[i]=m.classes[indmin(mydist)]
    end

    predict
    
end

function project(myLDA::PyObject,response::Array{Float64,2},stimulus::Array{Float64,1})

    predict=zeros(Float64,size(stimulus,1))

    classes=myLDA[:classes_]
    mydist=zeros(Float64,length(classes))
    mycentroids=myLDA[:transform](myLDA[:means_])
    
    for i=1:size(stimulus,1)
        xnew=myLDA[:transform](response[i,:])
        for j=1:length(classes)
            mydist[j]=norm(xnew-mycentroids[j,:])
        end
        predict[i]=classes[indmin(mydist)]
    end

    predict
    
end

function validate{C<:LDA, V<:Training}(m::decoder{C,V},response::Array{Float64,2},stimulus::Array{Float64,1})
    

end

function validate{C<:LDA, V<:LeaveOne}(m::decoder{C,V},response::Array{Float64,2},stimulus::Array{Float64,1})
    

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

function conmatrix{C<:classifier,V<:validation}(m::decoder{C,V},predict::Array{Float64,1},stimulus::Array{Float64,1})

    inds=collect(1:length(m.classes))

    mydict=Dict{Float64,Int64}([m.classes[i]=>inds[i] for i=1:length(m.classes)])

    #y is true, x is predicted
    #sum across columns should equal 1
    conmat=zeros(Float64,length(m.classes),length(m.classes))

    for i=1:length(stimulus)
        if length(find(stimulus[i].==m.classes))==0 #might be a stimulus in a val set that isn't in training
        else
            yind=mydict[stimulus[i]]
            xind=mydict[predict[i]]
            conmat[yind,xind]+=1
        end 
    end

    totals_real=sum(conmat,2)
    totals_predict=sum(conmat,1)
    
    for i=1:length(m.classes)
        conmat[i,:]=conmat[i,:]./totals_real[i]
    end

    p_real=squeeze(totals_real./sum(totals_real),2)
    p_predict=squeeze((totals_predict./sum(totals_predict))',2)
    
    (conmat, p_real, p_predict)
    
end
