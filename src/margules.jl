module margules

include("general.jl")

using Unitful
using NLsolve

∑(x)=sum(x)
γ(T,x,Λ)=(exp((Λ[1,2]+2*(Λ[2,1]-Λ[1,2])*x[1])*(x[2]^2)),exp((Λ[2,1]+2*(Λ[1,2]-Λ[2,1])*x[2])*(x[1]^2)))
P(T,x,CA,Λ,i) = x[i]*γ(T,x,Λ)[i]*general.P_Antoine(CA,T,i)
Pₜₒₜ(T,x,CA,Λ,n) = ∑(P(T,x,CA,Λ,i) for i=1:n)
y_calc(T,x,CA,Λ,i,n)=P(T,x,CA,Λ,i)/Pₜₒₜ(T,x,CA,Λ,n)

"""
    pxy(n,T,CA,Λ,xx;uni=u"Torr")

**Los campos de entrada son:**
- `n :: Int` Es el numero de componentes
- `T :: Float` Es la temperatura del sistema
- `CA :: n×3 Array{Float,2}` Constantes de Antoine para los `n` componentes
- `Λ :: n×n Array{Float,2}` Parámetros de interacción entre los componentes en grados kelvin
- `xx :: StepRange ó Vector` Son los puntos en el liquido donde se buscara el equilibrio con el vapor
- `uni` Son las unidades en las que se desea el resultado, por defecto esta en Torr

**Salida:**

Regresa las fracciones de liquido `x`,
las fracciones del vapor `y` y la presión `P` en Torr (`m` es la longitud de `xx`).

- `x :: n×m Array{Float64,2}` liquido
- `y :: n×m Array{Float64,2}` vapor
- `P :: m-element Array{Float64,1}` presión
"""
function pxy(n::Int,T,CA::Array,Λ::Array,xx;uni=u"Torr")
    T = T |> u"K"
    xx = general.nc(xx...)
    m = size(xx,2)
    y = Array{Float64}(undef,n,m)
    P = Vector{Float64}(undef,m)*uni
    for i=1:m
        x = xx[:,i]
        for j=1:n
            y[j,i] = y_calc(T,x,CA,Λ,j,n)
        end
        P[i] = Pₜₒₜ(T,x,CA,Λ,n) |> uni
    end
    
    return xx,y,P
end

"""
    txy(n,Pobj,CA,Λ,xx;uni=u"K")

**Los campos de entrada son:**
- `n :: Int` Es el numero de componentes
- `Pobj :: Float` Es la presión del sistema
- `CA :: n×3 Array{Float,2}` Constantes de Antoine para los `n` componentes
- `Λ :: n×n Array{Float,2}` Parámetros de interacción entre los componentes en grados kelvin
- `xx :: StepRange ó Vector` Son los puntos en el liquido donde se buscara el equilibrio con el vapor
- `uni` Son las unidades en las que se desea el resultado, por defecto esta en grados kelvin

**Salida:**

Regresa las fracciones de liquido `x`,
las fracciones del vapor `y` y la temperatura `T` en grados kelvin (`m` es la longitud de `xx`).

- `x :: n×m Array{Float64,2}` liquido
- `y :: n×m Array{Float64,2}` vapor
- `T :: m-element Array{Float64,1}` temperatura
"""
function txy(n::Int,Pobj,CA::Array,Λ::Array,xx;uni=u"K")
    Pobj = Pobj |> u"Torr"
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
        P2ₜₒₜ(T) = Pₜₒₜ(T,x,CA,Λ,n) - ustrip(u"Torr",Pobj)
        #println(i)
        Pₛₒₗ = nlsolve(n_ary(P2ₜₒₜ),[Tᵢ])
        Tᵢ = Pₛₒₗ.zero[1]
        T[i] = Tᵢ*u"K"
        for j = 1:n
            y[j,i] = y_calc(T[i],x,CA,Λ,j,n)
        end
    end
    
    return xx,y,T .|> uni
end

end  # module
