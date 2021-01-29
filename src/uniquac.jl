module uniquac

include("general.jl")
using Unitful
using NLsolve

∑(x)=sum(x)

τ(T::Quantity,Δu,i::Int,j::Int) = exp(-Δu[i,j]/(1.98721u"cal/(mol*K)"*T))
τ(T::Float64,Δu,i::Int,j::Int) = exp(-ustrip(u"cal/mol",Δu[i,j])/(1.98721*T))
V(x,r,i::Int,n::Int) = r[i]/∑(r[j]*x[j] for j=1:n)
F(x,q,i::Int,n::Int) = q[i]/∑(q[j]*x[j] for j=1:n)
γc(x,q,r,i::Int,n::Int) = 1-V(x,r,i,n)+log(V(x,r,i,n))-5*q[i]*(1-V(x,r,i,n)/F(x,q,i,n) +log(V(x,r,i,n)/F(x,q,i,n)))
γR(T,x,Δu,q,r,i::Int,n::Int) = q[i]*(1-log(∑(q[j]*x[j]*τ(T,Δu,j,i) for j=1:n)/∑(q[j]*x[j] for j=1:n))
-∑(q[j]*x[j]*τ(T,Δu,i,j)/∑(q[k]*x[k]*τ(T,Δu,k,j) for k=1:n) for j=1:n))
γ(T,x,Δu,q,r,i::Int,n::Int) = exp(γc(x,q,r,i,n) + γR(T,x,Δu,q,r,i,n))
P(T,x,CA,Δu,q,r,i::Int,n::Int) = x[i]*γ(T,x,Δu,q,r,i,n)*general.P_Antoine(CA,T,i)
Pₜₒₜ(T,x,CA,Δu,q,r,n::Int) = ∑(P(T,x,CA,Δu,q,r,i,n) for i=1:n)
y_calc(T,x,CA,Δu,q,r,i::Int,n::Int)=P(T,x,CA,Δu,q,r,i,n)/Pₜₒₜ(T,x,CA,Δu,q,r,n)

"""
    pxy(n,T,CA,Δu,q,r,xx;uni=u"Torr")

**Los campos de entrada son:**
- `n :: Int` Es el numero de componentes
- `T :: Float` Es la temperatura del sistema
- `CA :: n×3 Array{Float,2}` Constantes de Antoine para los `n` componentes
- `Δu :: n×n Array{Float,2}` Parametros de interacción entre los componentes en grados kelvin
- `q :: n-element Array{Float64,1}` área de Van der Waals
- `r :: n-element Array{Float64,1}` volumen de Van der Waals
- `xx :: StepRange ó Vector` Son los puntos en el liquido donde se buscara el equilibrio con el vapor
- `uni` Son las unidades en las que se desea el resultado, por defecto esta en Torr

**Salida:**

Regresa las fracciones de liquido `x`,
las fracciones del vapor `y` y la presión `P` en Torr (`m` es la longitud de `xx`).

- `x :: n×m Array{Float64,2}` liquido
- `y :: n×m Array{Float64,2}` vapor
- `P :: m-element Array{Float64,1}` presión
"""
function pxy(n::Int,T,CA::Array,Δu::Array,q::Vector,r::Vector,xx;uni=u"Torr")
    T = uconvert(u"K",T)
    xx = general.nc(xx...)
    m = size(xx,2)
    y = Array{Float64}(undef,n,m)
    P = Vector{Float64}(undef,m)*u"Torr"
    for i=1:m
        x = xx[:,i]
        for j=1:n
            y[j,i] = y_calc(T,x,CA,Δu,q,r,j,n)
        end
        P[i] = Pₜₒₜ(T,x,CA,Δu,q,r,n)
    end
    if uni==u"Torr"
        return xx,y,P
    else
        return xx,y,uconvert.(uni,P)
    end
end

"""
    txy(n,Pobj,CA,Δu,q,r,xx;uni=u"K")

**Los campos de entrada son:**
- `n :: Int` Es el numero de componentes
- `Pobj :: Float` Es la presión del sistema
- `CA :: n×3 Array{Float,2}` Constantes de Antoine para los `n` componentes
- `Δu :: n×n Array{Float,2}` Parametros de interacción entre los componentes en grados kelvin
- `q :: n-element Array{Float64,1}` área de Van der Waals
- `r :: n-element Array{Float64,1}` volumen de Van der Waals
- `xx :: StepRange ó Vector` Son los puntos en el liquido donde se buscara el equilibrio con el vapor
- `uni` Son las unidades en las que se desea el resultado, por defecto esta en grados kelvin

**Salida:**

Regresa las fracciones de liquido `x`,
las fracciones del vapor `y` y la temperatura `T` en grados kelvin (`m` es la longitud de `xx`).

- `x :: n×m Array{Float64,2}` liquido
- `y :: n×m Array{Float64,2}` vapor
- `T :: m-element Array{Float64,1}` temperatura
"""
function txy(n::Int,Pobj,CA::Array,Δu::Array,q::Vector,r::Vector,xx;uni=u"K")
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
        P2ₜₒₜ(T) = Pₜₒₜ(T,x,CA,Δu,q,r,n) - ustrip(u"Torr",Pobj)
        #println(i)
        Pₛₒₗ = nlsolve(n_ary(P2ₜₒₜ),[Tᵢ])
        Tᵢ = Pₛₒₗ.zero[1]
        T[i] = Tᵢ*u"K"
        for j = 1:n
            y[j,i] = y_calc(T[i],x,CA,Δu,q,r,j,n)
        end
    end
    if uni==u"K"
        return xx,y,T
    else
        return xx,y,uconvert.(uni,T)
    end
end
end  # module
