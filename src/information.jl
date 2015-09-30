


#=
Bias correction techniques

=#


#=
Quadratic extrapolation
Strong et al 1998
=#

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

#=
Decoding-based information

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
