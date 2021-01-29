module wilsonmod

include("general.jl")
using Unitful
using NLsolve

∑(x)=sum(x)
λ(T::Quantity,a,b,c,i::Int,j::Int) = a[i,j]u"cal/mol"+ b[i,j]u"cal/(mol*K)"*T + c[i,j]u"cal/(mol*K^2)"*T^2
λ(T::Float64,a,b,c,i::Int,j::Int) = a[i,j]+ b[i,j]*T + c[i,j]*T^2
Λ(T::Quantity,a,b,c,v,i::Int,j::Int) = (v[j]/v[i])*exp(-λ(T,a,b,c,i,j)/(1.98721u"cal/(mol*K)"*T))
Λ(T::Float64,a,b,c,v,i::Int,j::Int) = (v[j]/v[i])*exp(-λ(T,a,b,c,i,j)/(1.98721*T))
γ(T,x,a,b,c,v,i,n)=exp(-log(∑(x[j]*Λ(T,a,b,c,v,i,j) for j=1:n)) + 1 -∑((x[k]*Λ(T,a,b,c,v,k,i))/(∑(x[j]*Λ(T,a,b,c,v,k,j) for j=1:n)) for k=1:n))
P(T,x,CA,a,b,c,v,i,n) = x[i]*γ(T,x,a,b,c,v,i,n)*general.P_Antoine(CA,T,i)
Pₜₒₜ(T,x,CA,a,b,c,v,n) = ∑(P(T,x,CA,a,b,c,v,i,n) for i=1:n)
y_calc(T,x,CA,a,b,c,v,i,n) = P(T,x,CA,a,b,c,v,i,n)/Pₜₒₜ(T,x,CA,a,b,c,v,n)

"""
    pxy(n,T,CA,a,b,c,v,xx;uni=u"Torr")

**Los campos de entrada son:**
- `n :: Int` Es el numero de componentes
- `T :: Float` Es la temperatura del sistema
- `CA :: n×3 Array{Float,2}` Constantes de Antoine para los `n` componentes
- `a :: n×n Array{Float,2}` Parametros de interacción a en ``\\frac{cal}{mol}``
- `b :: n×n Array{Float,2}` Parametros de interacción b en ``\\frac{cal}{mol K}``
- `c :: n×n Array{Float,2}` Parametros de interacción c en ``\\frac{cal}{mol K^2}``
- `v :: n-element Array{Float,1}` Volumen molar
- `xx :: StepRange ó Vector` Son los puntos en el liquido donde se buscara el equilibrio con el vapor
- `uni` Son las unidades en las que se desea el resultado, por defecto esta en Torr

**Salida:**

Regresa las fracciones de liquido `x`,
las fracciones del vapor `y` y la presión `P` en Torr (`m` es la longitud de `xx`).

- `x :: n×m Array{Float64,2}` liquido
- `y :: n×m Array{Float64,2}` vapor
- `P :: m-element Array{Float64,1}` presión
"""
function pxy(n::Int,T,CA::Array,a::Array,b::Array,c::Array,v::Vector,xx;uni=u"Torr")
    T = uconvert(u"K",T)
    xx = general.nc(xx...)
    println("ddd",xx)
    m = size(xx,2)
    y = Array{Float64}(undef,n,m)
    P = Vector{Float64}(undef,m)*u"Torr"
    for i=1:m
        x = xx[:,i]
        for j=1:n
            y[j,i] = y_calc(T,x,CA,a,b,c,v,j,n)
        end
        P[i] = Pₜₒₜ(T,x,CA,a,b,c,v,n)
    end
    if uni==u"Torr"
        return xx,y,P
    else
        return xx,y,uconvert.(uni,P)
    end
end

"""
    txy(n,Pobj,CA,a,b,c,v,xx;uni=u"K")

**Los campos de entrada son:**
- `n :: Int` Es el numero de componentes
- `Pobj :: Float` Es la presión del sistema
- `CA :: n×3 Array{Float,2}` Constantes de Antoine para los `n` componentes
- `a :: n×n Array{Float,2}` Parametros de interacción a en ``\\frac{cal}{mol}``
- `b :: n×n Array{Float,2}` Parametros de interacción b en ``\\frac{cal}{mol K}``
- `c :: n×n Array{Float,2}` Parametros de interacción c en ``\\frac{cal}{mol K^2}``
- `v :: n-element Array{Float,1}` Volumen molar
- `xx :: StepRange ó Vector` Son los puntos en el liquido donde se buscara el equilibrio con el vapor
- `uni` Son las unidades en las que se desea el resultado, por defecto esta en grados kelvin

**Salida:**

Regresa las fracciones de liquido `x`,
las fracciones del vapor `y` y la temperatura `T` en grados kelvin (`m` es la longitud de `xx`).

- `x :: n×m Array{Float64,2}` liquido
- `y :: n×m Array{Float64,2}` vapor
- `T :: m-element Array{Float64,1}` temperatura
"""
function txy(n::Int,Pobj,CA::Array,a::Array,b::Array,c::Array,v::Vector,xx;uni=u"K")
    Pobj = uconvert(u"Torr",Pobj)
    if n==3
        xx = general.nc_ord(xx...)
    else
        xx = general.nc(xx...)
    end
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
        P2ₜₒₜ(T) = Pₜₒₜ(T,x,CA,a,b,c,v,n) - ustrip(u"Torr",Pobj)
        #println(i)
        Pₛₒₗ = nlsolve(n_ary(P2ₜₒₜ),[Tᵢ])
        Tᵢ = Pₛₒₗ.zero[1]
        T[i] = Tᵢ*u"K"
        for j = 1:n
            y[j,i] = y_calc(T[i],x,CA,a,b,c,v,j,n)
        end
    end
    if uni==u"K"
        return xx,y,T
    else
        return xx,y,uconvert.(uni,T)
    end
end
end  # module
