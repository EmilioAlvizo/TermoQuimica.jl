
module nrtl

include("general.jl")
using Unitful
using NLsolve

∑(x)=sum(x)

τ(T::Quantity,Δg,i::Int,j::Int) = Δg[i,j]/(1.98721u"cal/(mol*K)"*T)
τ(T::Float64,Δg,i::Int,j::Int) = ustrip(u"cal/mol",Δg[i,j])/(1.98721*T)
G(T,Δg,α,i::Int,j::Int) = exp(-α[i,j]*τ(T,Δg,i,j))
γ(T,x,Δg,α,i::Int,n::Int) = exp(∑(x[j]*G(T,Δg,α,j,i)*τ(T,Δg,j,i) for j=1:n)/∑(x[k]*G(T,Δg,α,k,i) for k=1:n)+ ∑((x[j]*G(T,Δg,α,i,j)/∑(x[k]*G(T,Δg,α,k,j) for k=1:n))*(τ(T,Δg,i,j)-∑(x[m]*G(T,Δg,α,m,j)*τ(T,Δg,m,j) for m=1:n)/∑(x[k]*G(T,Δg,α,k,j) for k=1:n)) for j=1:n))
P(T,x,CA,Δg,α,i::Int,n::Int) = x[i]*γ(T,x,Δg,α,i,n)*general.P_Antoine(CA,T,i)
Pₜₒₜ(T,x,CA,Δg,α,n::Int) = ∑(P(T,x,CA,Δg,α,i,n) for i=1:n)
y_calc(T,x,CA,Δg,α,i::Int,n::Int)=P(T,x,CA,Δg,α,i,n)/Pₜₒₜ(T,x,CA,Δg,α,n)

"""
    pxy(n,T,CA,Δg,α,xx;uni=u"Torr")

**Los campos de entrada son:**
- `n :: Int` Es el numero de componentes
- `T :: Float` Es la temperatura del sistema
- `CA :: n×3 Array{Float,2}` Constantes de Antoine para los `n` componentes
- `Δg :: n×n Array{Float,2}` Parametros de interacción entre los componentes en grados kelvin
- `α :: n×n Array{Float,2}` Parametro de no aleatoriedad
- `xx :: StepRange ó Vector` Son los puntos en el liquido donde se buscara el equilibrio con el vapor
- `uni` Son las unidades en las que se desea el resultado, por defecto esta en Torr

**Salida:**

Regresa las fracciones de liquido `x`,
las fracciones del vapor `y` y la presión `P` en Torr (`m` es la longitud de `xx`).

- `x :: n×m Array{Float64,2}` liquido
- `y :: n×m Array{Float64,2}` vapor
- `P :: m-element Array{Float64,1}` presión
"""
function pxy(n::Int,T,CA::Array,Δg::Array,α,xx;uni=u"Torr")
    T = uconvert(u"K",T)
    xx = general.nc(xx...)
    println("ddd",xx)
    m = size(xx,2)
    y = Array{Float64}(undef,n,m)
    P = Vector{Float64}(undef,m)*u"Torr"
    for i=1:m
        x = xx[:,i]
        for j=1:n
            y[j,i] = y_calc(T,x,CA,Δg,α,j,n)
        end
        P[i] = Pₜₒₜ(T,x,CA,Δg,α,n)
    end
    if uni==u"Torr"
        return xx,y,P
    else
        return xx,y,uconvert.(uni,P)
    end
end

"""
    txy(n,Pobj,CA,Δg,α,xx;uni=u"K")

**Los campos de entrada son:**
- `n :: Int` Es el numero de componentes
- `Pobj :: Float` Es la presión del sistema
- `CA :: n×3 Array{Float,2}` Constantes de Antoine para los `n` componentes
- `Δg :: n×n Array{Float,2}` Parametros de interacción entre los componentes en grados kelvin
- `α :: n×n Array{Float,2}` Parametro de no aleatoriedad
- `xx :: StepRange ó Vector` Son los puntos en el liquido donde se buscara el equilibrio con el vapor
- `uni` Son las unidades en las que se desea el resultado, por defecto esta en grados kelvin

**Salida:**

Regresa las fracciones de liquido `x`,
las fracciones del vapor `y` y la temperatura `T` en grados kelvin (`m` es la longitud de `xx`).

- `x :: n×m Array{Float64,2}` liquido
- `y :: n×m Array{Float64,2}` vapor
- `T :: m-element Array{Float64,1}` temperatura
"""
function txy(n::Int,Pobj,CA::Array,Δg::Array,α,xx;uni=u"K")
    Pobj = uconvert(u"Torr",Pobj)
    xx = general.nc(xx...)
    m = size(xx,2)
    y = Array{Float64}(undef,n,m)
    T = Vector{Float64}(undef,m)*u"K"
    Tᵢ = 300.0
    #temperatura inicial
    for i=1:n
        if xx[i,1]==1
            T[1] = general.T_Antoine(CA,Pobj,i)
            Tᵢ  = ustrip(u"K",T[1])
        end
    end
    for i = 1:m
        x = xx[:,i]
        P2ₜₒₜ(T) = Pₜₒₜ(T,x,CA,Δg,α,n) - ustrip(u"Torr",Pobj)
        #println(i)
        Pₛₒₗ = nlsolve(n_ary(P2ₜₒₜ),[Tᵢ])
        Tᵢ = Pₛₒₗ.zero[1]
        T[i] = Tᵢ*u"K"
        for j = 1:n
            y[j,i] = y_calc(T[i],x,CA,Δg,α,j,n)
        end
    end
    if uni==u"K"
        return xx,y,T
    else
        return xx,y,uconvert.(uni,T)
    end
end
end