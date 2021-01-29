using PolynomialRoots
using Unitful

∑(x)=sum(x)

const R = 8.3144621u"MPa*cm^3/(K*mol)"

function Van_Der_Waals_M(T,P,Tc,Pc,y)
    n = length(y)
    T = uconvert(u"K",T)
    P = uconvert(u"MPa",P)
    Tc = uconvert.(u"K",Tc)
    Pc = uconvert.(u"MPa",Pc)

    a = (27/64)*R^2 *Tc.^2 ./Pc
    b = R*Tc./(8*Pc)
    aa(i,j) =√(a[i]*a[j])
    aₘ = ∑(∑(y[i]*y[j]*aa(i,j) for j=1:n) for i=1:n)
    bₘ = ∑(y[i]*b[i] for i=1:n)

    A = aₘ*P/(R*T)^2
    B = bₘ*P/(R*T)
    

    𝒰 = -(1+B)
    𝒮 = A
    𝒯 = -A*B

    z1 = roots([𝒯, 𝒮, 𝒰, 1])

    index=findall(isreal,z1)
    if length(index)==3
        z = real(z1)
    end
    if length(index)==1
        z = real(z1[index[1]])
    end
    if length(index)==0
        for i=1:3
            if abs(imag(z1[i]))<=1e-5
                z = real(z1[i])
            end
        end
    end

    v = z*R*T/P

    U_rt = -aₘ./(R*T*v)
    H_rt = U_rt+z.-1
    G_rt = z.-1-log.(1 .-bₘ./v)-log.(z)-aₘ./(R*T*v)
    S_r = H_rt-G_rt
    A_rt = U_rt-S_r
    φ = exp.(G_rt)
    f = φ*P

    𝕌 = U_rt*R*T
    ℍ = H_rt*R*T
    𝕊 = S_r*R
    𝔾 = G_rt*R*T
    𝔸 = A_rt*R*T

    𝕌 = uconvert.(u"J/mol",𝕌)
    ℍ = uconvert.(u"J/mol",ℍ)
    𝕊 = uconvert.(u"J/(K*mol)",𝕊)
    𝔾 = uconvert.(u"J/mol",𝔾)
    𝔸 = uconvert.(u"J/mol",𝔸)

    ℍ = ustrip(ℍ)
    𝕌 = ustrip(𝕌)
    𝕊 = ustrip(𝕊)
    𝔾 = ustrip(𝔾)
    𝔸 = ustrip(𝔸)
    f = ustrip(f)
    v = ustrip(v)

    ℍ = round.(ℍ; digits=6)
    𝕌 = round.(𝕌; digits=6)
    𝕊 = round.(𝕊; digits=6)
    𝔾 = round.(𝔾; digits=6)
    𝔸 = round.(𝔸; digits=6)
    φ = round.(φ; digits=6)
    f = round.(f; digits=6)
    v = round.(v; digits=6)
    z = round.(z; digits=6)

    return ℍ,𝕌,𝕊,𝔾,𝔸,φ,f,z,v
end

function Redlich_Kwong_M(T,P,Tc,Pc,y)
    n = length(y)
    T = uconvert(u"K",T)
    P = uconvert(u"MPa",P)
    Tc = uconvert.(u"K",Tc)
    Pc = uconvert.(u"MPa",Pc)

    Tᵣ = T./Tc
    a = (0.42748*R^2 *Tc.^2.5)./Pc
    b = 0.08664*R*Tc./Pc
    
    aa(i,j) =√(a[i]*a[j])
    #aa(i,j) =(1-k)*√(a[i]*a[j])
    aₘ = ∑(∑(y[i]*y[j]*aa(i,j) for j=1:n) for i=1:n)
    bₘ = ∑(y[i]*b[i] for i=1:n)

    A = aₘ*P/(R^2 *T^2.5)
    B = bₘ*P/(R*T)
    𝒰 = -1
    𝒮 = A-B-B^2
    𝒯 = -A*B

    z1 = roots([𝒯, 𝒮, 𝒰, 1])

    index=findall(isreal,z1)
    if length(index)==3
        z = real(z1)
    end
    if length(index)==1
        z = real(z1[index[1]])
    end
    if length(index)==0
        for i=1:3
            if abs(imag(z1[i]))<=1e-5
                z = real(z1[i])
            end
        end
    end

    v = z*R*T/P

    U_rt = 3*aₘ*log.(z./(B.+z))/(2*bₘ*R*T^1.5)
    H_rt = U_rt+z.-1
    G_rt = z.-1-log.(z.-B)+(aₘ*log.(z./(z.+B)))./(bₘ*R*T^1.5)
    S_r = H_rt-G_rt
    A_rt = U_rt-S_r
    φ = exp.(G_rt)
    f = φ*P
    𝕌 = U_rt*R*T
    ℍ = H_rt*R*T
    𝕊 = S_r*R
    𝔾 = G_rt*R*T
    𝔸 = A_rt*R*T

    𝕌 = uconvert.(u"J/mol",𝕌)
    ℍ = uconvert.(u"J/mol",ℍ)
    𝕊 = uconvert.(u"J/(K*mol)",𝕊)
    𝔾 = uconvert.(u"J/mol",𝔾)
    𝔸 = uconvert.(u"J/mol",𝔸)

    ℍ = ustrip(ℍ)
    𝕌 = ustrip(𝕌)
    𝕊 = ustrip(𝕊)
    𝔾 = ustrip(𝔾)
    𝔸 = ustrip(𝔸)
    f = ustrip(f)
    v = ustrip(v)

    ℍ = round.(ℍ; digits=6)
    𝕌 = round.(𝕌; digits=6)
    𝕊 = round.(𝕊; digits=6)
    𝔾 = round.(𝔾; digits=6)
    𝔸 = round.(𝔸; digits=6)
    φ = round.(φ; digits=6)
    f = round.(f; digits=6)
    v = round.(v; digits=6)
    z = round.(z; digits=6)
    #return ℍ,𝕌,𝕊,z,v,φ,f
    return ℍ,𝕌,𝕊,𝔾,𝔸,φ,f,z,v

end # function

function Soave_Redlich_Kwong_M(T,P,Tc,Pc,y,𝜔)
    n = length(y)
    T = uconvert(u"K",T)
    P = uconvert(u"MPa",P)
    Tc = uconvert.(u"K",Tc)
    Pc = uconvert.(u"MPa",Pc)

    Tᵣ = T./Tc
    κ = 0.48508 .+ 1.55171*𝜔 -0.15613*𝜔.^2
    α = (1 .+ κ.*(1 .-Tᵣ.^.5)).^2
    ac = 0.42747*((R^2)*(Tc.^2)./Pc)

    a = ac.*α
    b = 0.08664*R*Tc./Pc
    A = a*P/(R*T)^2
    B = b*P/(R*T)
    
    aa(i,j) =√(a[i]*a[j])
    AA(i,j) =√(A[i]*A[j])
    #aa(i,j) =(1-k)*√(a[i]*a[j])
    aₘ = ∑(∑(y[i]*y[j]*aa(i,j) for j=1:n) for i=1:n)
    bₘ = ∑(y[i]*b[i] for i=1:n)
    Aₘ = aₘ*P/(R*T)^2
    Bₘ = bₘ*P/(R*T)

    𝒰 = -1
    𝒮 = Aₘ-Bₘ-Bₘ^2
    𝒯 = -Aₘ*Bₘ

    z1 = roots([𝒯, 𝒮, 𝒰, 1])
    

    index=findall(isreal,z1)
    if length(index)==3
        z = real(z1)
    end
    if length(index)==1
        z = real(z1[index[1]])
    end
    if length(index)==0
        for i=1:3
            if abs(imag(z1[i]))<=1e-5
                z = real(z1[i])
            end
        end
    end
    D = (-1/sqrt(T))*∑(∑(y[i]*y[j]*aa(i,j)*(κ[j]/sqrt(Tc[j]*α[j])+κ[i]/sqrt(Tc[i]*α[i]))/2 for j=1:n) for i=1:n)
    
    v = z*R*T/P

    U_rt = log(ustrip(R*T))*(aₘ-T*D)/(bₘ*R*T) .-log.(ustrip(R*T.+Bₘ*R*T./z))*(aₘ-T*D)/(bₘ*R*T)
    H_rt = z.-1 + U_rt
    A_rt = -log.(z) +(-R*T*log.(1 .-Bₘ./z)-aₘ*log.(1 .+Bₘ./z)./bₘ)/(R*T)
    S_r = -A_rt +U_rt 
    G_rt = H_rt-S_r
    
    φ = exp.(G_rt)
    f = φ*P
    𝕌 = U_rt*R*T
    ℍ = H_rt*R*T
    𝕊 = S_r*R
    𝔾 = G_rt*R*T
    𝔸 = A_rt*R*T

    𝕌 = uconvert.(u"J/mol",𝕌)
    ℍ = uconvert.(u"J/mol",ℍ)
    𝕊 = uconvert.(u"J/(K*mol)",𝕊)
    𝔾 = uconvert.(u"J/mol",𝔾)
    𝔸 = uconvert.(u"J/mol",𝔸)

    ℍ = ustrip(ℍ)
    𝕌 = ustrip(𝕌)
    𝕊 = ustrip(𝕊)
    𝔾 = ustrip(𝔾)
    𝔸 = ustrip(𝔸)
    f = ustrip(f)
    v = ustrip(v)

    ℍ = round.(ℍ; digits=6)
    𝕌 = round.(𝕌; digits=6)
    𝕊 = round.(𝕊; digits=6)
    𝔾 = round.(𝔾; digits=6)
    𝔸 = round.(𝔸; digits=6)
    φ = round.(φ; digits=6)
    f = round.(f; digits=6)
    v = round.(v; digits=6)
    z = round.(z; digits=6)
    #return ℍ,𝕌,𝕊,z,v,φ,f
    return ℍ,𝕌,𝕊,𝔾,𝔸,φ,f,z,v

end # function

function Peng_Robinson_M(T,P,Tc,Pc,y,𝜔)
    n = length(y)
    T = uconvert(u"K",T)
    P = uconvert(u"MPa",P)
    Tc = uconvert.(u"K",Tc)
    Pc = uconvert.(u"MPa",Pc)

    Tᵣ = T./Tc
    κ = 0.37464 .+ 1.54226*𝜔 -0.26993*𝜔.^2
    α = (1 .+ κ.*(1 .-Tᵣ.^.5)).^2
    ac = 0.457235529*((R^2)*(Tc.^2)./Pc)

    a = ac.*α
    b = 0.077796074*R*Tc./Pc
    A = a*P/(R*T)^2
    B = b*P/(R*T)


    aa(i,j) =√(a[i]*a[j])
    AA(i,j) =√(A[i]*A[j])
    #aa(i,j) =(1-k)*√(a[i]*a[j])
    aₘ = ∑(∑(y[i]*y[j]*aa(i,j) for j=1:n) for i=1:n)
    bₘ = ∑(y[i]*b[i] for i=1:n)
    Aₘ = aₘ*P/(R*T)^2
    Bₘ = bₘ*P/(R*T)

    𝒰 = +Bₘ-1
    𝒮 = Aₘ-2*Bₘ-3*Bₘ^2
    𝒯 = -(Aₘ*Bₘ-Bₘ^2-Bₘ^3)

    z = roots([𝒯, 𝒮, 𝒰, 1])
    println(z)

    index=findall(x->abs(imag(x))>1e-5,z)
    deleteat!(z,index)
    z = real(z)
    index=findall(x->x<0,z)
    deleteat!(z,index)
    if length(z)==1
        z=z[1]
    end
    println(z)
    
    v = z*R*T/P

    D = (-1/√T)*∑(∑(y[i]*y[j]*aa(i,j)*(κ[j]/√(Tc[j]*α[j])+κ[i]/√(Tc[i]*α[i]))/2 for j=1:n) for i=1:n)
    H_rt = -1 .+z+(acoth(√2).+atanh.((-1 .+Bₘ./z)/(√2)))*(-aₘ+T*D)/(√2*bₘ*R*T)
    U_rt = H_rt-(z.-1)
    A_rt = -log.(z.-Bₘ)-aₘ*(acoth(√2).+atanh.((-1 .+Bₘ./z)/(√2)))/(√2*bₘ*R*T)
    S_r = U_rt-A_rt
    G_rt = H_rt-S_r
  
    φ = exp.(G_rt)
    f = φ*P
    𝕌 = U_rt*R*T
    ℍ = H_rt*R*T
    𝕊 = S_r*R
    𝔾 = G_rt*R*T
    𝔸 = A_rt*R*T

    𝕌 = uconvert.(u"J/mol",𝕌)
    ℍ = uconvert.(u"J/mol",ℍ)
    𝕊 = uconvert.(u"J/(K*mol)",𝕊)
    𝔾 = uconvert.(u"J/mol",𝔾)
    𝔸 = uconvert.(u"J/mol",𝔸)

    ℍ = ustrip(ℍ)
    𝕌 = ustrip(𝕌)
    𝕊 = ustrip(𝕊)
    𝔾 = ustrip(𝔾)
    𝔸 = ustrip(𝔸)
    f = ustrip(f)
    v = ustrip(v)

    ℍ = round.(ℍ; digits=6)
    𝕌 = round.(𝕌; digits=6)
    𝕊 = round.(𝕊; digits=6)
    𝔾 = round.(𝔾; digits=6)
    𝔸 = round.(𝔸; digits=6)
    φ = round.(φ; digits=6)
    f = round.(f; digits=6)
    v = round.(v; digits=6)
    z = round.(z; digits=6)
    return ℍ,𝕌,𝕊,𝔾,𝔸,φ,f,z,v
end 

Peng_Robinson_M(100u"K",10u"atm",[155.5,388.75]u"K",[10.1,25.25]u"atm",[.5,.5],[.155,.3875])


Redlich_Kwong_M(100u"°F",278u"psi",[190.6,305.4,369.8]u"K",[45.4,48.2,41.9]u"atm",[.1,.67,.53],1)

Redlich_Kwong_M(30u"°C",25.5u"atm",[304.2,364.9]u"K",[72.9,45.45]u"atm",[.5,.5])
Van_Der_Waals_M(30u"°C",25.5u"atm",[304.2,364.9]u"K",[72.9,45.45]u"atm",[.5,.5])


Redlich_Kwong_M(55u"°C",279.005u"atm",[126.2,154.6]u"K",[33.5,49.8]u"atm",[.79,.21])
Van_Der_Waals_M(55u"°C",279.005u"atm",[126.2,154.6]u"K",[33.5,49.8]u"atm",[.79,.21])
Soave_Redlich_Kwong_M(55u"°C",279.005u"atm",[126.2,154.6]u"K",[33.5,49.8]u"atm",[.79,.21],[0.04,0.021])


Redlich_Kwong_M(400u"°F",30u"atm",[305.1,369.8,425.2]u"K",[48.2,41.9,37.5]u"atm",[1/3,1/3,1/3])
Van_Der_Waals_M(400u"°F",30u"atm",[305.1,369.8,425.2]u"K",[48.2,41.9,37.5]u"atm",[1/3,1/3,1/3])
Soave_Redlich_Kwong_M(400u"°F",30u"atm",[305.1,369.8,425.2]u"K",[48.2,41.9,37.5]u"atm",[1/3,1/3,1/3],[.098,.152,.198])
Peng_Robinson(400u"°F",30u"atm",[305.1,369.8,425.2]u"K",[48.2,41.9,37.5]u"atm",[1/3,1/3,1/3],[.098,.152,.198])

Redlich_Kwong_M(150u"K",2.8u"MPa",[126.1,190.6,425.2,540.3,591.8]u"K",[3.394,4.604,3.797,2.736,4.109]u"MPa",[.600056,.399944,6.81e-12,1.42e-14,3.41e-14])
Van_Der_Waals_M(150u"K",2.8u"MPa",[126.1,190.6,425.2,540.3,591.8]u"K",[3.394,4.604,3.797,2.736,4.109]u"MPa",[.600056,.399944,6.81e-12,1.42e-14,3.41e-14])
Soave_Redlich_Kwong_M(150u"K",2.8u"MPa",[126.1,190.6,425.2,540.3,591.8]u"K",[3.394,4.604,3.797,2.736,4.109]u"MPa",[.600056,.399944,6.81e-12,1.42e-14,3.41e-14],[.04,.011,.193,.349,.264])
Peng_Robinson_M(400u"K",2.8u"MPa",[126.1,190.6,425.2,540.3,591.8]u"K",[3.394,4.604,3.797,2.736,4.109]u"MPa",[.600056,.399944,6.81e-12,1.42e-14,3.41e-14],[.04,.011,.193,.349,.264])


Soave_Redlich_Kwong_M(200u"K",10u"Torr",[155.5,388.75,570.1666666]u"K",[10.1,25.25,37.033333333]u"atm",[1/3,1/3,1/3],[.155,.3875,.568333333])


ff=[0.121729801 0.1435876	0.41687813	0.711599775	0.582443308;
0.1435876	0.162856597	0.469570996	0.805827602	0.649337699;
0.41687813	0.469570996	1.352245643	2.322817807	1.866413123;
0.711599775	0.805827602	2.322817807	3.987043948	3.210696736;
0.582443308	0.649337699	1.866413123	3.210696736	2.568751717]

ss=[0.600055682,0.399944318,6.80916e-12,1.42249e-14,3.41449e-14]
[0.245930877 0.48954508 4.362334064 12.17868034 9.220793434]*[0.245930877,0.48954508,4.362334064,12.17868034,9.220793434]
[0.245930877,0.48954508,4.362334064,12.17868034,9.220793434]*[0.245930877 0.48954508 4.362334064 12.17868034 9.220793434]


(ff*ss)
(ff*ss)'

(ff*ss)'*ss
(ff*ss).*ss