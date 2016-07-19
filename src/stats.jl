
export anova, anova2, anova_repeated, circ_anova

function anova(x,alpha=.05)
    
    df_b=size(x,2)-1
    df_w=length(x)-size(x,2)
    df_t=df_b+df_w
       
    myd=FDist(df_b,df_w)
    
    means=mean(x,1)
    
    grand_mean=mean(x)
    
    ss_total=0.0
    for i=1:length(x)
        ss_total+=((x[i]-grand_mean)*(x[i]-grand_mean))
    end
    
    ss_w=0.0
    for i=1:size(x,2)
        for j=1:size(x,1)
            ss_w += ((x[j,i]-means[i])*(x[j,i]-means[i]))
        end
    end

    ss_b=ss_total-ss_w
    
    us_b=ss_b/df_b    
    us_w=ss_w/df_w
    
    F=us_b/us_w
    
    p=ccdf(myd,F)
    sig=false
    if p<alpha
        sig=true
    end

    p
end

function anova2(x,a,b)
    
    group_a=unique(a)
    
    group_b=unique(b)
    
    m_a=zeros(Int64,length(unique(a)))
    m_b=zeros(Int64,length(unique(b)))
    
    df_a=length(m_a)-1
    df_b=length(m_b)-1
    
    df_ab=df_a*df_b
    
    df_w=length(x)-length(m_a)*length(m_b)
    
    df_t=length(x)-1
    
    sumx=zeros(Float64,length(m_a),length(m_b))
    sumx2=zeros(Float64,size(sumx))
    
    n=zeros(Int64,size(sumx))
    
    SS=zeros(Float64,size(sumx))
    
    sumx2n=zeros(Float64,size(sumx))
    a_ind=1
    b_ind=1
    for i=1:length(x)
        for j=1:length(group_a)
            if group_a[j]==a[i]
                a_ind=group_a[j]
                break
            end
        end
        for j=1:length(group_b)
            if group_b[j]==b[i]
                b_ind=group_b[j]
                    break
            end
        end
        
        sumx[a_ind,b_ind]+=x[i]
        sumx2[a_ind,b_ind]+=(x[i]*x[i])
        
        n[a_ind,b_ind]+=1
    end
    
    for i=1:length(sumx)
       
        sumx2n[i]=(sumx[i]*sumx[i])/n[i]
        
        SS[i]=sumx2[i]-sumx2n[i]
        
    end
    
    SS_w=sum(SS)
    
    SS_b=sum([sum(sumx[:,i])^2/sum(n[:,i]) for i=1:size(sumx,2)])-(sum(sumx))^2/sum(n)
    
    SS_a=sum([sum(sumx[i,:])^2/sum(n[i,:]) for i=1:size(sumx,1)])-(sum(sumx))^2/sum(n)

    SS_ab=sum(sumx2n)-sum([sum(sumx[:,i])^2/sum(n[:,i]) for i=1:size(sumx,2)])-sum([sum(sumx[i,:])^2/sum(n[i,:]) for i=1:size(sumx,1)])+(sum(sumx))^2/sum(n)
    
    SS_t=SS_w+SS_a+SS_b+SS_ab
    
    MS_a=SS_a/df_a
    
    MS_b=SS_b/df_b
    
    MS_ab=SS_ab/df_ab
    
    MS_w=SS_w/df_w
    
    F_a=MS_a/MS_w
    
    F_b=MS_b/MS_w
    
    F_ab=MS_ab/MS_w
    
    myda=FDist(df_a,df_w)
    pa=ccdf(myda,F_a)
    
    mydb=FDist(df_b,df_w)
    pb=ccdf(mydb,F_b)
    
    mydab=FDist(df_ab,df_w)
    pab=ccdf(mydab,F_ab)
    
    (pa,pb,pab)
end

function anova_repeated(x)
    
    means=mean(x,1) 
    
    g_mean=mean(x)
    
    ss_b=0.0
    
    for i=1:size(x,2)
        ss_b+=(means[i]-g_mean)^2
    end
    
    ss_b*=size(x,1)
    
    ss_w=0.0
    
    for i=1:size(x,2)
        for j=1:size(x,1)
            ss_w+=(x[j,i]-means[i])^2
        end
    end
    
    means_s=mean(x,2)
    ss_s=0.0
    
    for i=1:size(x,1)
        ss_s+=(means_s[i]-g_mean)^2
    end
    
    ss_s*=size(x,2)
    
    ss_e=ss_w-ss_s
    
    us_b=ss_b/(size(x,2)-1)
    
    us_e=ss_e/((size(x,1)-1)*(size(x,2)-1))
    
    F=us_b/us_e
    
    myd=FDist(size(x,2-1),size(x,1)-1)
    
    
    ccdf(myd,F)
end

function circ_anova(x,alpha=.05)
    
    A=zeros(Complex,size(x))
    
    for i=1:length(x)
        A[i] = exp(im*x[i])
    end
    
    R = abs(mean(A))    
    N=length(A)
    mu = angle(mean(A))
    
    #Group means
    Ri=abs(mean(A,1))

    #Group concentrations
    ki=zeros(Float64,size(A,2))
    for i = 1:size(A,2)
        ki[i] = Concentration(A[:,i],Ri[i])
    end
    
    SSw=N
    for i=1:size(A,2)
        SSw -= size(A,1)*Ri[i]
    end
    SSb = 0.0
    for i=1:size(A,2)
        SSb += size(A,1)*Ri[i]
    end
    SSb -= N*R

    F = (N-size(A,2))/(size(A,2)-1)*SSb/SSw
    k=sum([size(A,1)*ki[i] for i=1:size(A,2)])/N
    if 2 < k && k < 10
        #Improve X² approximation for moderate concentration (Stephens, 1969)
        F = F * (1+3/(8*k))
    end
    
    myd=FDist(size(A,2)-1,N-size(A,2))
    
    p=ccdf(myd,F)
    sig=false
    if p<alpha
        sig=true
    end

    (sig,p)
end

function Concentration(A,r_bar)
    
    n=length(A)
    
    if r_bar < 0.53
        kappa = 2*r_bar+r_bar^3+5*r_bar^5/6
    elseif r_bar < 0.85
        kappa = -0.4+1.39*r_bar+0.43/(1-r_bar)
    else
        kappa = 1/(r_bar^3-4*r_bar^2+3*r_bar)
    end

    #Correction for small samples
    if n <= 15
        if kappa < 2
            kappa = max([kappa-2/(n*kappa); 0.0])
        else
            kappa = (n-1)^3*kappa/(n^3+n)
        end
    end
    kappa
end
