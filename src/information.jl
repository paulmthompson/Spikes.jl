
export MI_bias_QE

#=
Bias correction techniques

=#


#=
Quadratic extrapolation
Strong et al 1998

I=I_raw + a/n_s + b/n_s^2

Compute MI with number of repetitions (n_s) equal to total, total/2, and total/4
Plot result against 1/n_s, and fit quadratic polynomial by least squares
Intercept is true information
=#

function MI_bias_QE{C<:classifier,V<:validation}(m::decoder{C,V},response::Array{Float64,2},stimulus::Array{Float64,1})

    N_s=length(realval)
    
    myinds=1:N_s

    samplepoints=[N_s; N_s/2; N_s/4]
    
    MI=zeros(Float64,3)   

    predict=project(m,response,stimulus)
    (conmat, p_s, p_r)=conmatrix(m,predict,stimulus)

    MI[1]=mutual_information(conmat,p_s,p_r)

    for i=2:3

        theseinds=sample(myinds, samplepoints[i], replace=false)

        stimulus_2=stimulus[theseinds]
        response_2=response[theseinds,:]
        predict_2=project(m,response_2,stimulus_2)
        
        (conmat2, p_s2, p_r2)=conmatrix(m,predict_2,stimulus_2)

        MI[i]=mutual_information(conmat2,p_s2,p_r2)

    end

    I_corrected=polyfit(1./samplepoints, MI, 2)

end

#=
Panzeri-Treves (PT)
Panzeri and Treves 1996
=#

#=
Nemenman-Shafee-Bialek (NSB)
Nemenman et al 2002
=#

#=
James-Stein shrinkage
Hausser and Strimmer 2008
=#

function mutual_information(p_rs::Array{Float64,2},p_s::Array{Float64,1},p_r::Array{Float64,1})

    MI=0.0
    
    for i=1:size(p_rs,1)
        for j=1:size(p_rs,2)
            MI+=p_rs[i,j]*p_s[i]*log2(p_rs[i,j]/(p_r[j]))      
        end
    end

    MI
    
end

function polyfit(x::Array{Float64,1}, y::Array{Float64,1}, n::Int64)
  A = [ x[i]^p for i = 1:length(x), p = 0:n ]
  A \ y
end
