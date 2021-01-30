module general
using Unitful

"""
    P_Antoine(CA, T, i)

```math
\\log_{10}{P(mmHg)}=A-\\frac{B}{C+T(°C)}
```

**Los campos de entrada son:**

- `CA :: n×3 Array{Float,2}` Constantes de Antoine para los `n` componentes
- `T :: Float` Es la temperatura del sistema en °C
- `i :: Int` Componente `i` al que se quiere conocer su presión

**Salida:**

Regresa la presión `P` en mmHg

- `P :: m-element Array{Float64,1}` presión
"""
function P_Antoine(CA::Array, T::Unitful.Temperature, i::Int)::Unitful.Pressure
    A,B,C = CA[i, :]
    T = T |> u"K"
    P = u"Torr"*10^(A-B/((T/1u"K" -273.15)+C))
    return P
end

"""
    T_Antoine(CA, P, i)

```math
T(°C)=-C-\\frac{B}{log_{10}{P(mmHg)}-A}
```

**Los campos de entrada son:**

- `CA :: n×3 Array{Float,2}` Constantes de Antoine para los `n` componentes
- `P :: Float` Es la presión del sistema en mmHg
- `i :: Int` Componente `i` al que se quiere conocer su temperatura

**Salida:**

Regresa la temperatura `T` en °C

- `T :: m-element Array{Float64,1}` temperatura
"""
function T_Antoine(CA::Array, P::Unitful.Pressure, i::Int)::Unitful.Temperature
    A,B,C = CA[i, :]
    P = P |> u"Torr"
    (-B/(log10(P/1u"Torr")-A) -C +273.15)*u"K"
end


function P_Antoine(CA::Array, T, i::Int)::Float64
    A,B,C = CA[i, :]
    P = 10^(A - B / ((T - 273.15) + C))
    return P
end


function T_Antoine(CA::Array, P, i::Int)::Float64
    A,B,C = CA[i, :]
    (-B / (log10(P) - A) - C + 273.15)
end


function nelemen(n::Int, ñ::Int)::Int
    m = Array{Int}(undef,ñ,n)
    for i=1:n, j=1:ñ
        if i == 1 || j == 1
            m[j, i] = 1
        else
            m[j, i] = m[j, i - 1] + m[j - 1, i]
        end
    end
    return m[ñ, n]
end
function nc(xx::StepRangeLen...)
    kk = collect(Iterators.product(xx...))
    sk = sum.(kk)
    m = length(xx)
    ñ = length(xx[1])
    n = length(sk)
    nn = nelemen(m + 1,ñ)
    x = Array{Float64}(undef,m + 1,nn)
    k = 1
    for i=1:n
        if sk[i] <= 1
            for j=1:m
                x[j, k] = kk[i][j]
            end
            x[m + 1, k] = 1 - sk[i]
            k +=1
        end
    end
    return x
end
function nc(xx::Array...)
    m = length(xx)
    ñ = length(xx[1])
    x = Array{Float64}(undef,m + 1,ñ)
    for i=1:ñ
        sk = 0
        for j=1:m
            x[j, i] = xx[j][i]
            sk += xx[j][i]
        end
        if sk <= 1
            x[m + 1, i] = 1 - sk
        else
            println("La suma de las primeras variables es mayor a 1")
        end
    end
    return x
end
function nc_ord(xx::StepRangeLen...)
    kk = collect(Iterators.product(xx...))
    sk = sum.(kk)
    m = length(xx)
    ñ = length(xx[1])
    n = length(sk)
    nn = nelemen(m + 1,ñ)
    x = Array{Float64}(undef,m + 1,nn)
    k = 1
    iter = 1
    for i=1:ñ
        if iter == 1
            for h=1:ñ
                if sk[i, h] <= 1
                    for j=1:m
                        x[j, k] = kk[i, h][j]
                    end
                    x[m + 1, k] = 1 - sk[i, h]
                    k +=1
                end
            end
            iter = 0
        else
            for h=reverse(1:ñ)
                if sk[i, h] <= 1
                    for j=1:m
                        x[j, k] = kk[i, h][j]
                    end
                    x[m + 1, k] = 1 - sk[i, h]
                    k +=1
                end
            end
            iter = 1
        end
    end
    return x
end

end  # module
