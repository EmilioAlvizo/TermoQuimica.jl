using Unitful
using BenchmarkTools
using NLsolve

push!(LOAD_PATH,"/Users/emili/AppData/Local/Programs/Julia/Julia-1.4.1/aaa/TermoQuimica/src")
using fricción
using tuberias

1+5

f_kw(x; y = 3) = x + y
function fi(f,s,a) 
    xx=y=f*s/a
    println(xx)
    println(y)
end
function fi(y) 
    xx=y=f*s/a
    println(xx)
    println(y)
end

fi(2,5,3)
fi(2.56)
@btime "Rennels"===4
@btime "Rennels"==4
@btime f_kw(2)
@btime f_kw(2, y = 3)
@btime f(2)
@btime f(2,5)

mutable struct flujo
    L::typeof(1.0u"m")
    sh::String
    Dᵢ::typeof(1.0u"mm")
    Dₒ::typeof(1.0u"mm")
    ΔD::typeof(1.0u"mm")
    NPS::Float64
    Q::typeof(1.0u"l/hr")
    ṁ::typeof(1.0u"kg/hr")
    v::typeof(1.0u"m/s")
    ρ::typeof(1.0u"kg/m^3")
    𝜈::typeof(1.0u"m^2/s")
    μ::typeof(1.0u"cP")
    Re::Float64
    ε::typeof(1.0u"mm")
    f::Float64
    hf::typeof(1.0u"m")
end

function flujo(L,sh;Dᵢ,Dₒ,ΔD,NPS,Q,ṁ,v,ρ,𝜈,μ,Re,ε,f,hf)
    flujo(L,sh,Dᵢ,Dₒ,ΔD,NPS,Q,ṁ,v,ρ,𝜈,μ,Re,f,hf)
end

function flujo(L,sh;ṁ,ρ,μ,Dₒ,ε)
    NPS,Dᵢ,Dₒ,ΔD = tuberia(sh,Do=Dₒ)
    Q = ṁ/ρ
    v = Q/(π/4*Dᵢ^2)
    𝜈 = μ/ρ
    Re = upreferred(v*Dᵢ/𝜈)
    f = factor_friccion(Re,ε/Dᵢ)
    hf = upreferred(f*(L/Dᵢ)*(v^2 /(2*9.80665u"m/s^2")))
    flujo(L,sh,Dᵢ,Dₒ,ΔD,NPS,Q,ṁ,v,ρ,𝜈,μ,Re,ε,f,hf)
end

function flujo(L,sh;Q,ρ,𝜈,Dᵢ,ε)
    NPS,Dᵢ,Dₒ,ΔD = tuberia(sh,Di=Dᵢ)
    ṁ = Q*ρ
    v = Q/(π/4*Dᵢ^2)
    μ = 𝜈*ρ
    Re = upreferred(v*Dᵢ/𝜈)
    f = factor_friccion(Re,ε/Dᵢ)
    hf = upreferred(f*(L/Dᵢ)*(v^2 /(2*9.80665u"m/s^2")))
    flujo(L,sh,Dᵢ,Dₒ,ΔD,NPS,Q,ṁ,v,ρ,𝜈,μ,Re,ε,f,hf)
end

flujo(500u"m","40",0.2u"m^3/s",900u"kg/m^3",0.00001u"m^2/s",Dᵢ=200u"mm",ε=.26u"mm")

flujo(120u"m","80",40000u"kg/hr",12.493325597325967u"kg/m^3",2.655716090653184e-05u"Pa*s",Dₒ=150u"mm",ε=.0018u"inch")


.0018u"inch"
flujo(5u"m",5u"m",5u"m",5u"m",5u"m",5u"m",5u"m",5u"m",5u"m",5u"m",5u"m",5u"m",5u"m",5u"m")
flujo(5u"cm",5u"m",5u"m",5u"m",5u"m",5u"m",5u"m",5u"m",5u"m",5u"m",5u"m",5u"m",5u"m",5u"m")

tuberia("40",Do=200u"mm")

factor_friccion(1, 1E-4)
factor_friccion(128000, .0013)
factor_friccion(1, 1E-4,función=Wood)
factor_friccion(128000, .0013,función=Moody)

laminar(100)

Reynolds(D,v,ρ,μ) = upreferred(ρ*v*D/μ)
Reynolds(D,v,𝜈) = upreferred(v*D/𝜈)

Reynolds(2.5u"m", 0.25u"m/s", 1.1613u"kg/m^3", 1.9E-5u"Pa*s")
Reynolds(2.5u"m", 0.25u"m/s", 1.1613u"kg/m^3", 1.9E-5u"cP")
Reynolds(2.5u"m", 0.25u"m/s", 1.1613u"kg/m^3", 1.9E-5u"N*s/m^2")
Reynolds(2.5u"m", 0.25u"m/s", 1.636e-05u"m^2/s")


function hf(L::Unitful.Length,D::Unitful.Length,v::Unitful.Velocity,ρ::Unitful.Density,μ::Unitful.DynamicViscosity,eD;método::Function=Colebrook)
    upreferred(factor_friccion(Reynolds(D,v,ρ,μ), eD,función=método)*(L/D)*(v^2 /(2*9.80665u"m/s^2")))
end
function hf(L::Unitful.Length,D::Unitful.Length,v::Unitful.Velocity,𝜈::Unitful.KinematicViscosity,eD;método::Function=Colebrook)
    upreferred(factor_friccion(Reynolds(D,v,𝜈), eD,función=método)*(L/D)*(v^2 /(2*9.80665u"m/s^2")))
end

hf(400u"m",150u"mm",0.0157u"m/s",1.1613u"kg/m^3", 1.9E-3u"Pa*s",0)
hf(500u"m",200u"mm",6.4u"m/s",.00001u"m^2/s",.0013)


Unitful.Velocity