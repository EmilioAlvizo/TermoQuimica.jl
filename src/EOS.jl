module EOS
using PolynomialRoots
using Unitful

∑(x)=sum(x)

const R = 8.31446261815324u"MPa*cm^3/(K*mol)"

@doc raw"""
    Van_Der_Waals(T,P,Tc,Pc)

**Ecuación de la forma:**

```math
P=\frac{RT}{ν-b}-\frac{a}{ν^2}
```
Donde:
```math
\begin{aligned}
a&=\frac{27(RT_c)^2}{64P_c} & b&=\frac{RT_c}{8P_c}\\
\end{aligned}
```
```math
\frac{U-U^{ig}}{RT}=\frac{a}{RTν}
```
```math
\frac{G-G^{ig}}{RT}=z-1-\frac{a}{RTν}-\log_{10} (z-B)
```
**Los campos de entrada son:**
- `T :: Float` Es la temperatura del sistema
- `P :: Float` Es la presión del sistema
- `Tc :: Float` Es la temperatura critica del sistema
- `Pc :: Float` Es la presión critica del sistema
- `𝜔 :: Float` Es el factor acéntrico 

**Salida:**

Regresa un conjunto de propiedades termodinámicas.

- `ℍ :: Float` Entalpía residual 
- `𝕌 :: Float` Energía interna residual
- `𝕊 :: Float` Entropía residual
- `𝔾 :: Float` Energía libre de Gibbs residual
- `𝔸 :: Float` Energía libre de helmholtz residual
- `φ :: Float` Coeficiente de fugacidad 
- `f :: Float` Fugacidad
- `z :: Float` Factor de compresibilidad 
- `ν :: Float` Volumen molar

"""
function Van_Der_Waals(T,P,Tc,Pc)
    T = uconvert(u"K",T)
    P = uconvert(u"MPa",P)
    Tc = uconvert(u"K",Tc)
    Pc = uconvert(u"MPa",Pc)

    a = (27/64)*R^2 *Tc^2 /Pc
    b = R*Tc/(8*Pc)
    A = a*P/(R*T)^2
    B = b*P/(R*T)

    𝒰 = -(1+B)
    𝒮 = A
    𝒯 = -A*B

    z = roots([𝒯, 𝒮, 𝒰, 1])
    index=findall(x->abs(imag(x))>1e-5,z)
    deleteat!(z,index)
    z = real(z)
    index=findall(x->x<0,z)
    deleteat!(z,index)
    if length(z)==1
        z=z[1]
    end

    v = z*R*T/P

    U_rt = -a./(R*T*v)
    H_rt = U_rt+z.-1
    G_rt = z.-1-log.(z.-B)-a./(R*T*v)
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

@doc raw"""
    Redlich_Kwong(T,P,Tc,Pc)

**Ecuación de la forma:**

```math
P=\frac{RT}{ν-b}-\frac{a}{(T^{1.5} ν(ν+b))}
```
Donde:
```math
\begin{aligned}
a&=\frac{0.42748 R^2 T_c^{2.5}}{P_c} & b&=\frac{0.08664 RT_c}{P_c}\\
\end{aligned}
```
```math
\frac{U-U^{ig}}{RT}=\frac{3a\log_{10} (1+z/B)}{2bR T^{1.5}}
```
```math
\frac{G-G^{ig}}{RT}=z-1+\frac{a\log_{10} (1+z/B)}{bR T^{1.5}}-\log_{10} (z-B)
```

**Los campos de entrada son:**
- `T :: Float` Es la temperatura del sistema
- `P :: Float` Es la presión del sistema
- `Tc :: Float` Es la temperatura critica del sistema
- `Pc :: Float` Es la presión critica del sistema
- `𝜔 :: Float` Es el factor acéntrico 

**Salida:**

Regresa un conjunto de propiedades termodinámicas.

- `ℍ :: Float` Entalpía residual 
- `𝕌 :: Float` Energía interna residual
- `𝕊 :: Float` Entropía residual
- `𝔾 :: Float` Energía libre de Gibbs residual
- `𝔸 :: Float` Energía libre de helmholtz residual
- `φ :: Float` Coeficiente de fugacidad 
- `f :: Float` Fugacidad
- `z :: Float` Factor de compresibilidad 
- `v :: Float` Volumen molar
"""
function Redlich_Kwong(T,P,Tc,Pc)
    T = uconvert(u"K",T)
    P = uconvert(u"MPa",P)
    Tc = uconvert(u"K",Tc)
    Pc = uconvert(u"MPa",Pc)

    Tᵣ = T/Tc
    a = (0.42748*R^2 *Tc^2.5)/Pc
    b = 0.08664*R*Tc/Pc

    A = a*P/(R^2 *T^2.5)
    B = b*P/(R*T)
    𝒰 = -1
    𝒮 = A-B-B^2
    𝒯 = -A*B

    z = roots([𝒯, 𝒮, 𝒰, 1])
    index=findall(x->abs(imag(x))>1e-5,z)
    deleteat!(z,index)
    z = real(z)
    index=findall(x->x<0,z)
    deleteat!(z,index)
    if length(z)==1
        z=z[1]
    end

    v = z*R*T/P
    U_rt = 3*a*log.(z./(B.+z))/(2*b*R*T^1.5)
    H_rt = U_rt+z.-1
    G_rt = z.-1-log.(z.-B)+(a*log.(z./(z.+B)))./(b*R*T^1.5)
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

@doc raw"""
    Soave_Redlich_Kwong(T,P,Tc,Pc,𝜔)

**Ecuación de la forma:**

```math
P=\frac{RT}{V-b}-\frac{aα}{V(V+b)}
```
Donde:
```math
\begin{aligned}
a&=\frac{0.42748 R^2 T_c^{2}}{P_c} & b&=\frac{0.08664 RT_c}{P_c} & α&=(1 + κ(1 -T_r^{0.5}))^2\\
\end{aligned}
```
```math
κ=0.48+1.574w-0.176w^2
```
```math
\frac{U-U^{ig}}{RT}=-\frac{(aα+a\sqrt{T_rα}κ)\log_{10} (1+B/z)}{bRT}
```
```math
\frac{S-S^{ig}}{R}=\log_{10}(z-B)-\frac{a\sqrt{T_rα}κ\log_{10}(1+B/z)}{bRT}
```

**Los campos de entrada son:**
- `T :: Float` Es la temperatura del sistema
- `P :: Float` Es la presión del sistema
- `Tc :: Float` Es la temperatura critica del sistema
- `Pc :: Float` Es la presión critica del sistema
- `𝜔 :: Float` Es el factor acéntrico 

**Salida:**

Regresa un conjunto de propiedades termodinámicas.

- `ℍ :: Float` Entalpía residual 
- `𝕌 :: Float` Energía interna residual
- `𝕊 :: Float` Entropía residual
- `𝔾 :: Float` Energía libre de Gibbs residual
- `𝔸 :: Float` Energía libre de helmholtz residual
- `φ :: Float` Coeficiente de fugacidad 
- `f :: Float` Fugacidad
- `z :: Float` Factor de compresibilidad 
- `v :: Float` Volumen molar
"""
function Soave_Redlich_Kwong(T,P,Tc,Pc,𝜔)
    T = uconvert(u"K",T)
    P = uconvert(u"MPa",P)
    Tc = uconvert(u"K",Tc)
    Pc = uconvert(u"MPa",Pc)

    Tᵣ = T/Tc
    κ = 0.48 + 1.574*𝜔 -0.176*𝜔^2
    α = (1 + κ*(1 -Tᵣ^.5))^2
    ac = 0.42748*((R^2)*(Tc^2)/Pc)
    a = ac*α
    b = 0.08664*R*Tc/Pc

    A = a*P/(R*T)^2
    B = b*P/(R*T)
    𝒰 = -1
    𝒮 = A-B-B^2
    𝒯 = -A*B

    z = roots([𝒯, 𝒮, 𝒰, 1])
    index=findall(x->abs(imag(x))>1e-5,z)
    deleteat!(z,index)
    z = real(z)
    index=findall(x->x<0,z)
    deleteat!(z,index)
    if length(z)==1
        z=z[1]
    end

    v = z*R*T/P
    U_rt = -((a.+ac*√Tᵣ*√α*κ)*log.(1 .+B./z))/(b*R*T)
    H_rt = U_rt+z.-1
    S_r = log.(z.-B)-(ac*√Tᵣ*√α*κ*log.(1 .+B./z))/(b*R*T)
    G_rt = H_rt-S_r
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

@doc raw"""
    Peng_Robinson(T,P,Tc,Pc,𝜔)

**Ecuación de la forma:**

```math
P=\frac{RT}{V-b}-\frac{aα}{V^2 +2bV+b^2}
```
Donde:
```math
\begin{aligned}
a&=\frac{0.45724 R^2 T_c^{2}}{P_c} & b&=\frac{0.0778 RT_c}{P_c} & α&=(1 + κ(1 -T_r^{0.5}))^2\\
\end{aligned}
```
```math
κ=0.37464 + 1.54226*w -0.26992*w^2
```
```math
\frac{U-U^{ig}}{RT}=-(aα+a\sqrt{T_rα}κ)*(\coth^{-1}{\sqrt{2}}+\tanh^{-1}\frac{(-1+B/z)/\sqrt{2}}{bRT\sqrt{2}})
```
```math
\frac{G-G^{ig}}{R}=z-1-\log(z-B)-aα*(\coth^{-1}{\sqrt{2}}+\tanh^{-1}\frac{(B-z)/z\sqrt{2}}{bRT\sqrt{2}})
```

**Los campos de entrada son:**
- `T :: Float` Es la temperatura del sistema
- `P :: Float` Es la presión del sistema
- `Tc :: Float` Es la temperatura critica del sistema
- `Pc :: Float` Es la presión critica del sistema
- `𝜔 :: Float` Es el factor acéntrico 

**Salida:**

Regresa un conjunto de propiedades termodinámicas.

- `ℍ :: Float` Entalpía residual 
- `𝕌 :: Float` Energía interna residual
- `𝕊 :: Float` Entropía residual
- `𝔾 :: Float` Energía libre de Gibbs residual
- `𝔸 :: Float` Energía libre de helmholtz residual
- `φ :: Float` Coeficiente de fugacidad 
- `f :: Float` Fugacidad
- `z :: Float` Factor de compresibilidad 
- `v :: Float` Volumen molar
"""
function Peng_Robinson(T,P,Tc,Pc,𝜔)
    T = uconvert(u"K",T)
    P = uconvert(u"MPa",P)
    Tc = uconvert(u"K",Tc)
    Pc = uconvert(u"MPa",Pc)

    Tᵣ = T/Tc
    κ = 0.37464 + 1.54226*𝜔 -0.26992*𝜔^2
    α = (1 + κ*(1 -Tᵣ^.5))^2
    ac = 0.45724*((R^2)*(Tc^2)/Pc)
    a = ac*α
    b = 0.0778*R*Tc/Pc

    𝒰 = b*P/(R*T) -1
    𝒮 = a*P/(R*T)^2 -3*(b*P/(R*T))^2 -2*b*P/(R*T)
    𝒯 = (b*P/(R*T))^3 +(b*P/(R*T))^2 -a*P*b*P/(R*T)^3

    z = roots([𝒯, 𝒮, 𝒰, 1])
    index=findall(x->abs(imag(x))>1e-5,z)
    deleteat!(z,index)
    z = real(z)
    index=findall(x->x<0,z)
    deleteat!(z,index)
    if length(z)==1
        z=z[1]
    end

    v = z*R*T/P
    A1 = a*P/(R*T)^2
    B = b*P/(R*T)

    U_rt = -(a.+ac*√Tᵣ*√α*κ)*(acoth(√2).+atanh.((-1 .+B./z)./√2))/(√2*b*R*T)
    H_rt = U_rt+z.-1
    G_rt = z.-1-log.(z.-B)-a*(acoth(√2).+atanh.((B.-z)./(z*√2)))/(√2*b*R*T)
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

    z = roots([𝒯, 𝒮, 𝒰, 1])
    index=findall(x->abs(imag(x))>1e-5,z)
    deleteat!(z,index)
    z = real(z)
    index=findall(x->x<0,z)
    deleteat!(z,index)
    if length(z)==1
        z=z[1]
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

    z = roots([𝒯, 𝒮, 𝒰, 1])
    index=findall(x->abs(imag(x))>1e-5,z)
    deleteat!(z,index)
    z = real(z)
    index=findall(x->x<0,z)
    deleteat!(z,index)
    if length(z)==1
        z=z[1]
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

    z = roots([𝒯, 𝒮, 𝒰, 1])
    index=findall(x->abs(imag(x))>1e-5,z)
    deleteat!(z,index)
    z = real(z)
    index=findall(x->x<0,z)
    deleteat!(z,index)
    if length(z)==1
        z=z[1]
    end

    D = (-1/sqrt(T))*∑(∑(y[i]*y[j]*aa(i,j)*(κ[j]/sqrt(Tc[j]*α[j])+κ[i]/sqrt(Tc[i]*α[i]))/2 for j=1:n) for i=1:n)
    
    v = z*R*T/P

    U_rt = log(ustrip(R*T))*(aₘ-T*D)/(bₘ*R*T) .-log.(ustrip(R*T.+Bₘ*R*T./z))*(aₘ-T*D)/(bₘ*R*T)
    H_rt = z.-1 + U_rt
    A_rt = -log.(z) +(-R*T*log.(1 .-Bₘ./z)-aₘ*log.(1 .+Bₘ./z)/bₘ)/(R*T)
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

"""
    Peng_Robinson_M(T,P,Tc,Pc,y,𝜔)

**Ecuación de la forma:**

```math
P=\\frac{RT}{V-b}-\\frac{aα}{V^2 +2bV+b^2}
```

**Los campos de entrada son:**
- `T :: Float` Es la temperatura del sistema
- `P :: Float` Es la presión del sistema
- `Tc :: Float` Es la temperatura critica del sistema
- `Pc :: Float` Es la presión critica del sistema
- `𝜔 :: Float` Es el factor acéntrico 

**Salida:**

Regresa un conjunto de propiedades termodinámicas.

- `ℍ :: Float` Entalpía residual 
- `𝕌 :: Float` Energía interna residual
- `𝕊 :: Float` Entropía residual
- `𝔾 :: Float` Energía libre de Gibbs residual
- `𝔸 :: Float` Energía libre de helmholtz residual
- `φ :: Float` Coeficiente de fugacidad 
- `f :: Float` Fugacidad
- `z :: Float` Factor de compresibilidad 
- `v :: Float` Volumen molar
"""
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
    index=findall(x->abs(imag(x))>1e-5,z)
    deleteat!(z,index)
    z = real(z)
    index=findall(x->x<0,z)
    deleteat!(z,index)
    if length(z)==1
        z=z[1]
    end

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


end

