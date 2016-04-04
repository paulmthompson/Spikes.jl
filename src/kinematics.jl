
export get_rt_d

function get_rt_d(kin::Array{Kinematic,1},thres::Float64)

    rt=zeros(Int64,length(kin))
    
    for i=1:length(kin)
        c1=(mean(kin[i].px[1:10]),mean(kin[i].py[1:10]))

        t_rt=0
        for j=10:length(kin[i].px)
            if (kin[i].px[j]-c1[1])^2 + (kin[i].py[j]-c1[2])^2 > thres*thres
                t_rt=j
                break
            end
        end

        rt[i]=t_rt
    end
    rt
end

function plot_kinx(kin::Array{Kinematic,1},tp::Array{Int64,1},myind::Array{Int64,1},l::Int64,ax)
    y=get_kinx(kin,tp,myind)
    (myci,mymean)=plot_kin(y,l)
    shaded_line(mymean,myci,ax)
end

function plot_kinxt(kin::Array{Kinematic,1},tp::Array{Int64,1},myind::Array{Int64,1},l::Int64,ax)

    x=[collect(1:length(tp[i]-9:length(kin[i].vx))) for i=1:length(kin[myind])]
    y=get_kinx(kin,tp,myind)
    xy=[hcat(x[i],y[i]) for i=1:length(x)]

    c1=C.LineCollection(xy,color="black",linewidth=1.0,alpha=.2)
    ax[:add_collection](c1)
    ax[:set_xlim]([0.0, l])
    ax[:set_ylim]([-1.0, 1.0])

    nothing
end

function plot_kiny(kin::Array{Kinematic,1},tp::Array{Int64,1},myind::Array{Int64,1},l::Int64,ax)
    y=get_kiny(kin,tp,myind)
    (myci,mymean)=plot_kin(y,l)
    shaded_line(mymean,myci,ax)
end

function plot_kinyt(kin::Array{Kinematic,1},tp::Array{Int64,1},myind::Array{Int64,1},l::Int64,ax)

    x=[collect(1:length(tp[i]-9:length(kin[i].vy))) for i=1:length(kin[myind])]
    y=get_kiny(kin,tp,myind)
    xy=[hcat(x[i],y[i]) for i=1:length(x)]

    c1=C.LineCollection(xy,color="black",linewidth=1.0,alpha=.2)
    ax[:add_collection](c1)
    ax[:set_xlim]([0.0, l])
    ax[:set_ylim]([-1.0, 1.0])

    nothing   
end

function plot_kin(y::Array{Array{Float64,1},1},l::Int64)

    myci=zeros(Float64,l,2)
    mymean=zeros(Float64,l)
    for i=1:l
        test=zeros(Float64,0)
        for j=1:length(y)
            if length(y[j])>i
                push!(test,y[j][i])
            end
        end
        #temp=ci(OneSampleTTest(test),tail = :both)
        temp=quantile(test,[.25,.75])
        myci[i,1]=temp[1]
        myci[i,2]=temp[2]
        mymean[i]=mean(test)
    end
    (myci,mymean)
end

function get_kinx(kin::Array{Kinematic,1},tp::Array{Int64,1},myind::Array{Int64,1})
    vx=[kin[i].vx[tp[i]-9:end] for i=1:length(kin[myind])]
end

function get_kiny(kin::Array{Kinematic,1},tp::Array{Int64,1},myind::Array{Int64,1})
    vy=[kin[i].vy[tp[i]-9:end] for i=1:length(kin[myind])]
end

function shaded_line(y::Array{Float64,1},myci::Array{Float64,2},ax)
    ax[:plot](collect(1:length(y)),y,color="blue")
    ax[:fill_between](collect(1:length(y)),myci[:,1],myci[:,2],color="blue",alpha=0.5)
    nothing
end
