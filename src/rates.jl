

export rate_session, ses_mean, ses_std, rate_event, rate_window, rate_window_pop, zscore, zscore_pop, psth_per_con,plot_raster, plot_psth, plot_psth_raster, plot_raster_psth, plot_pop_psth, plot_pop_psth_cb, psth_pop


#=
Temporal rate for entire experimental session
=#

rate_session(myrate::rate_bin,n::Int64)=rate_session_(myrate.spikes[n],myrate.binsize)

function rate_session_(spikes::SpikeTrain,binsize::Float64)
    hist(spikes.ts,collect(spikes.trials[1].time:binsize:spikes.trials[end].time))[2]
end

#=
Temporal rate within time window
=#

function spike_count_win(rate::SpikeTrain,t1::Float64,t2::Float64)

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

#=
Mean rate within window
=#

#=
Trial-averaged rate within window
=#

#=
Trial-averaged temporal rate within window
=#

rate_window(myrate::rate_bin,event_id::Int,ts,n::Int64)=rate_window(myrate.spikes[n],event_id,ts,myrate.binsize)

function rate_window(spikes::SpikeTrain,event_id::Int,ts,binsize::Float64)
    spikehist(myrate,ind,ts,n)./myrate.binsize
end

function spikehist(myrate::rate,ind::Int64,ts::FloatRange{Float64},n::Int64)
    hist(myrate.spikes[n].ts[myrate.spikes[n].trials[ind].inds]-myrate.spikes[n].center[ind,1],ts[1]:myrate.binsize:ts[end])[2]
end

function spikehist(myrate::rate,ind::Int64,ts::Array{Float64,1},n::Int64)
    hist(myrate.spikes[n].ts[myrate.spikes[n].trials[ind].inds]-myrate.spikes[n].center[ind,1],ts)[2]
end

function rate_event(myrate::rate_bin,inds::Array{Int64,1},ts::FloatRange{Float64},n::Int64)
    collect_ts(myrate,inds,ts,n)./(ts.step/ts.divisor)
end



function rate_event(myrate::rate_bin,inds::Array{Int64,1},ts::Array{Float64,1},n::Int64)
    collect_ts(myrate,inds,ts,n)./(ts[2]-ts[1])
end

function rate_event(myrate::rate_bin,inds::Array{Int64,1},ts::Array{Float64,1},n::Int64,output::Array{Float64,1})
    collect_ts(myrate,inds,ts,n,output)
    @inbounds for i=1:length(output)
        output[i]=output[i]/(ts[i+1]-ts[i])
    end
    nothing
end

function rate_event(myrate::rate_bin,inds::Array{Int64,1},ts::FloatRange{Float64},n::Int64,output::Array{Float64,1})
    collect_ts(myrate,inds,collect(ts),n,output)
end

#=
TODO
Shimazaki H. and Shinomoto S., A method for selecting the bin size of a time histogram
=#

#=
Kernel Density Methods

=#
function rate_session(myrate::rate_KD,n::Int64)
    binned=hist(myrate.spikes[n].ts,myrate.spikes[n].trials[1].time:myrate.binsize:myrate.spikes[n].trials[end].time)[2]
    imfilter_gaussian(binned,[myrate.kern/myrate.binsize])./myrate.binsize
end

function rate_window(myrate::rate_KD,ind::Int64,time::FloatRange{Float64},n::Int64)
    binned=spikehist(myrate,ind,time,n)
    imfilter_gaussian(binned,[myrate.kern/myrate.binsize])./myrate.binsize 
end

function rate_event(myrate::rate_KD,inds::Array{Int64,1},time::FloatRange{Float64},n::Int64)
    binned=collect_ts(myrate,inds,time,n)
    imfilter_gaussian(binned,[myrate.kern/myrate.binsize])./myrate.binsize   
end

#=
Shimazaki H. and Shinomoto S., Kernel Bandwidth Optimization in Spike Rate Estimation
=#
function rate_event_kdopt(spikes::SpikeTrain,inds::Array{Int64,1},time::Array{Float64,1})  
    spikes_temp=sort(vcat([spikes.ts[spikes.trials[i].inds]-spikes.center[i,1] for i in inds]...))

    isi=diff(spikes_temp)
    dt=minimum(isi)

    #t=linspace(time[1],time[2],round(Int64,min((time[2]-time[1])/dt+.5,1000)))
    t=linspace(time[1],time[2],1000)
    
    dt=minimum(diff(t))

    binned=hist(spikes_temp, [t-dt/2; t[end]+dt/2])[2]
    
    L=length(binned)
    N=length(spikes_temp)
    
    wmin=2*dt
    wmax=time[2]-time[1]
    imax=20
    
    w=zeros(Float64,imax)
    c=zeros(Float64,imax)
    optw=0.0
    tol=1e-5
    
    phi=.5*(sqrt(5)+1)
    
    a=ilogexp(wmin)
    b=ilogexp(wmax)
    
    c1= (phi-1) * a + (2-phi) * b
    c2= (2 - phi)*a*(phi-1) *b
    
    f1=cost_function(binned, N, logexp(c1), dt)
    f2=cost_function(binned, N, logexp(c1), dt)
    
    k=1
    
    while ((abs(b-a)) > (tol * (abs(c1) + abs(c2)))) & (k<imax)
        if f1<f2
            b=c2
            c2=c1
            c1 = (phi-1)*a + (2-phi) * b
            f2=f1
            f1=cost_function(binned,N,logexp(c1),dt)
            w[k]=logexp(c1)
            c[k]=f1
            optw=logexp(c1)
            
        else
            a=c1
            c1=c2
            c2=(2-phi)*a+(phi-1)*b
            f1=f2
            f2=cost_function(binned,N,logexp(c1),dt)
            w[k]=logexp(c2)
            c[k]=f2
            optw=logexp(c2)
        end
        k=k+1
    end
    
    yh=abs(imfilter_gaussian(binned, [optw/dt]))./(length(inds)*dt)

    (yh, optw/dt)
end

function ilogexp(x::Float64)
   if x < 1e2
        y = log(exp(x)-1)
    else
        y = x
    end 
    
    y
end

function logexp(x::Float64)
    if x < 1e2 
        y = log(1 + exp(x))
    else
        y = x
    end
    y
end

function cost_function(binned::Array{Int64,1},N::Int64,w::Float64,dt::Float64)
    #The cost function
    #Cn(w) = sum_{i,j} int k(x - x_i) k(x - x_j) dx - 2 sum_{i~=j} k(x_i - x_j) 
    
    yh=abs(imfilter_gaussian(binned, [w/dt]))
    
    C=sum(yh.*2)*dt - 2 * sum(yh.*binned) * dt + 2 / sqrt(2*pi) / w / N
    C = C * N * N
    
    C   
end

#=
Helper Methods
=#


#Mean firing rate for entire session
function ses_mean(myrate::rate,n::Int64)
    mean(rate_session(myrate,n))
end

#STD of firing rate for whole session
function ses_std(myrate::rate,n::Int64)
    std(rate_session(myrate,n))
end

#
function rate_trials(myrate::rate,inds::Array{Int64,1},time::FloatRange{Float64},n::Int64)

    spikes_temp=zeros(Float64,length(inds),length(time[1]:myrate.binsize:time[end])-1)
    for i=1:length(inds)
        spikes_temp[i,:]=rate_window(myrate,inds[i],time,n)    
    end
    spikes_temp   
end

function spikehist(myrate::rate,ind::Int64,ts::FloatRange{Float64},n::Int64)
    hist(myrate.spikes[n].ts[myrate.spikes[n].trials[ind].inds]-myrate.spikes[n].center[ind,1],ts[1]:myrate.binsize:ts[end])[2]
end

function spikehist(myrate::rate,ind::Int64,ts::Array{Float64,1},n::Int64)
    hist(myrate.spikes[n].ts[myrate.spikes[n].trials[ind].inds]-myrate.spikes[n].center[ind,1],ts)[2]
end

function collect_ts(myrate::rate,inds::Array{Int64,1},ts::FloatRange{Float64},n::Int64,output=zeros(Float64,length(ts)-1))
    collect_ts(myrate,inds,collect(ts),n,output)
end

function collect_ts(myrate::rate,inds::Array{Int64,1},ts::Array{Float64,1},n::Int64,output=zeros(Float64,length(ts)-1))
    @inbounds for i=1:length(output)
        output[i]=0.0
    end
   @inbounds for i in inds
        count=1
        mycent=myrate.spikes[n].center[i,1]
        for j=1:length(myrate.spikes[n].trials[i].inds)
            myt=myrate.spikes[n].ts[myrate.spikes[n].trials[i].inds[j]]-mycent
            if myt>ts[count]
                k=count+1
                while k < length(ts)+1
                    if myt<ts[k]
                        output[k-1]+=1.0
                        count=k-1
                        break
                    end
                    k+=1
                end
            end
        end
    end
    @inbounds for i=1:length(output)
        output[i]=output[i]/length(inds)
    end
    output
end

function zscore(myrate::rate,inds::Array{Int64,1},ts::FloatRange{Float64},n::Int64)
    mypsth=rate_event(myrate,inds,ts,n)
    wholetrain=rate_session(myrate,n)
    (mypsth-mean(wholetrain))./std(wholetrain)
end

function zscore(myrate::rate,inds::Array{Int64,1},ts::Array{Float64,1},n::Int64)
    mypsth=rate_event(myrate,inds,ts,n)
    wholetrain=rate_session(myrate,n)
    (mypsth-mean(wholetrain))./std(wholetrain)
end

function rate_window_pop(myrate::rate,ind::Int64,ts::FloatRange{Float64})
     rate_window_pop(myrate,ind,collect(ts))
end

function rate_window_pop(myrate::rate,ind::Int64,ts::Array{Float64,1})
    raster=zeros(Float64,length(time[1]:myrate.binsize:time[end])-1,length(myrate.spikes))

    for i=1:size(raster,2)
        raster[:,i]=rate_window(myrate,ind,time,i)
    end
    raster 
end

zscore_pop(myrate::rate,inds::Array{Int64,1},ts::FloatRange{Float64})=zscore_pop(myrate,inds,collect(ts))

function zscore_pop(myrate::rate,inds::Array{Int64,1},ts::Array{Float64,1})
    raster=zeros(Float64,length(myrate.spikes),length(ts)-1)

    for k=1:size(raster,1)
        raster[k,:]=zscore(myrate,inds,ts,k)
    end   
    raster
end

psth_pop(myrate::rate,inds::Array{Int64,1},ts::FloatRange{Float64},raster=zeros(Float64,length(myrate.spikes),length(ts)-1))=psth_pop(myrate,inds,collect(ts),raster)

function psth_pop(myrate::rate,inds::Array{Int64,1},ts::Array{Float64,1},raster=zeros(Float64,length(myrate.spikes),length(ts)-1))   
    tempmat=zeros(Float64,length(ts)-1)
    for k=1:size(raster,1)
        rate_event(myrate,inds,ts,k,tempmat)
        @inbounds for j=1:length(tempmat)
            raster[k,j]=tempmat[j]
        end
    end
    raster
end

function total_spikes(myrate::rate,inds::Array{Int64,1},n::Int64)
    mysize=0
    for i in inds
        @inbounds mysize+=length(myrate.spikes[n].trials[i].inds)
    end
    mysize
end

function total_spikes(myrate::rate,inds::Array{Int64,1},ts::Array{Float64,1},n::Int64)

    mysize=0
    @inbounds for i in inds
        mycent=myrate.spikes[n].center[i,1]
        for j in myrate.spikes[n].trials[i].inds
            if myrate.spikes[n].ts[j]-mycent > ts[end]
                break
            elseif myrate.spikes[n].ts[j]-mycent > ts[1]
                mysize+=1
            end
        end
    end
    mysize
end
    
function psth_per_con(myrate::rate,ts::Array{Float64,1},inds::Array{Int64,1},n::Int64,cons::Array{Int64,1})

    psths=zeros(Float64,length(ts)-1,length(cons))

    for i=1:length(cons)
        con_inds=find(inds.==cons[i])
        psths[:,i]=rate_event(myrate,con_inds,ts,n)
    end

    psths
end

#=
Plotting
=#

#=
Raster
=#

function plot_raster(myrate::rate,inds::Array{Int64,1},ts::Array{Float64,1},n::Int64,ax; event_ts=zeros(Float64,0))

    plot_raster(myrate,inds,ts,n,ax)

    if length(event_ts)>0

        xy=[zeros(Float64,2,2) for i=1:(size(event_ts,2)*length(inds))]
            
        mycount=1.0
        myinds=1
        for i in inds
            for j=1:size(event_ts,2)
                xy[myinds][1,1]=event_ts[i]-.015
                xy[myinds][1,2]=mycount+.5
                xy[myinds][2,1]=event_ts[i]+.015
                xy[myinds][2,2]=mycount+.5
                myinds+=1
            end
            mycount+=1.0
        end
        c2=C.LineCollection(xy,color="black",linewidth=2.0)
        ax[:add_collection](c2)
    end
    nothing
end

function plot_raster(myrate::rate,inds::Array{Int64,1},ts::Array{Float64,1},n::Int64,ax)

    mycount=1.0
    myspikes=1
    xy=[zeros(Float64,2,2) for i=1:total_spikes(myrate,inds,ts,n)]
    @inbounds for i in inds
        mycenter=myrate.spikes[n].center[i,1] 
        for j in myrate.spikes[n].trials[i].inds
            if myrate.spikes[n].ts[j]-mycenter > ts[end]
                break
            elseif myrate.spikes[n].ts[j] - mycenter > ts[1]
                xy[myspikes][1,1]=myrate.spikes[n].ts[j]-mycenter
                xy[myspikes][1,2]=mycount
                xy[myspikes][2,1]=myrate.spikes[n].ts[j]-mycenter
                xy[myspikes][2,2]=mycount+1.0
                myspikes+=1
            end
        end
        mycount+=1.0
    end
    c1=C.LineCollection(xy,color="black",linewidth=1.0)
    ax[:add_collection](c1)
    ax[:set_yticks]([])
    ax[:set_ylabel]("")
    ax[:set_ylim]([1,mycount])
    ax[:set_xticks]([])
    ax[:set_xticklabels]([])
    ax[:set_xlim]([ts[1],ts[end]])
    nothing
end

#=
PSTH
=#

plot_psth(myrate::rate,ts::FloatRange{Float64},inds::Array{Int64,1},n::Int64)=plot_psth(myrate,collect(ts),inds,n)

plot_psth(myrate::rate,ts::FloatRange{Float64},inds::Array{Int64,1},n::Int64,ax)=plot_psth(myrate,collect(ts),inds,n,ax)

function plot_psth(myrate::rate,ts::Array{Float64,1},inds::Array{Int64,1},n::Int64)

    (fig,ax)=subplots(1,1)

    plot_psth(myrate,ts,inds,n,ax)

    (fig,ax)
end

function plot_psth(myrate::rate,ts::Array{Float64,1},inds::Array{Int64,1},n::Int64,ax)

    psth=rate_event(myrate,inds,ts,n)

    ylimit=round(Int64,maximum(psth)/10)*10
    ax[:set_yticks]([0, ylimit])
    ax[:set_yticklabels]([0,ylimit],size=6)
    ax[:set_ylabel]("Rate (spikes/s)", size=6)

    ax[:set_xticks]([])
    ax[:set_xticklabels]([])
    ax[:bar](2:length(ts),psth,color="black")
    ax[:set_ylim]([0, ylimit])
    nothing
end

plot_psth_raster(myrate::rate,ts::FloatRange{Float64},inds::Array{Int64,1},n::Int64)=plot_psth_raster(myrate,collect(ts),inds,n)

function plot_psth_raster(myrate::rate,ts::Array{Float64,1},inds::Array{Int64,1},n::Int64)

    (fig,ax)=subplots(2,1)
    subplots_adjust(hspace=0.0)

    plot_psth(myrate,ts,inds,n,ax[1,1])
    plot_raster(myrate,inds,ts,n,ax[2,1])
    ax[2,1][:set_xlabel]("Time (s)", size=8)
    (fig,ax)
end

function plot_psth_raster(myrate::rate,ts::Array{Float64,1},inds::Array{Int64,1},n::Int64,ax_x::Int64,ax_y::Int64,ax;event_ts=zeros(Float64,0))
    plot_psth(myrate,ts,inds,n,ax[ax_x,ax_y])
    plot_raster(myrate,inds,ts,n,ax[ax_x+1,ax_y],event_ts=event_ts)
    nothing
end

function plot_raster_psth(myrate::rate,ts::FloatRange{Float64},inds::Array{Int64,1},n::Int64)

    (fig,ax)=subplots(2,1)
    subplots_adjust(hspace=0.0)
    
    plot_raster(myrate,inds,ts,n,ax[1,1])
    plot_psth(myrate,ts,inds,n,ax[2,1])
    ax[2,1][:set_xlabel]("Time (s)", size=8)
    (fig,ax)
end

function plot_pop_psth(myrate::rate,ts::FloatRange{Float64},inds::Array{Int64,1},lims::NTuple{2,Float64})

    (fig,ax)=subplots(1,1)
    myplot=plot_pop_psth(myrate,ts,inds,lims,ax)
    (fig,ax)
end

function plot_pop_psth(myrate::rate,ts::FloatRange{Float64},inds::Array{Int64,1},lims::NTuple{2,Float64},ax)

    popth=zscore_pop(myrate,inds,ts)

    ax[:set_yticks]([])
    ax[:set_yticklabels]([])
    ax[:set_ylabel]("Neuron", size=8)

    ax[:set_xlabel]("Time (s)", size=8)

    myplot=ax[:pcolor](popth,cmap="jet",vmin=lims[1],vmax=lims[2])
end

function plot_pop_psth_cb(myrate::rate,ts::FloatRange{Float64},inds::Array{Int64,1},lims::NTuple{2,Float64})
    (fig,ax)=subplots(1,1)
    myplot=plot_pop_psth(myrate,ts,inds,lims,ax)
    cbar_ax = fig[:add_axes]([.85,.15,.02,.7])
    fig[:colorbar](myplot,cax=cbar_ax,label="Z score",ticks=[lims[1],lims[2]])
    (fig,ax)
end
    
