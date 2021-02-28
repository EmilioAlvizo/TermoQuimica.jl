module saft

using Unitful
using NLsolve

include("general.jl")

∑(x)=sum(x)
const q01, q11, q21 = [3.7039,-0.3226,0.6907]u"Å"
const q02, q12, q22 = [0.06233,-0.02236,-0.01563]u"mol/g"
const q03, q13, q23 = [150.03,80.68,38.96]u"K"
const MCH4 = 16.043u"g/mol"
const R = 8.31446261815324u"MPa*cm^3/(K*mol)"
const Nₐ = 6.02214076e23u"1/(mol)"
const kᵦ = R/Nₐ

const a0 = [0.9105631445,0.6361281449,2.6861347891,-26.547362491,97.759208784,-159.59154087,91.297774084]
const a1 = [-0.3084016918,0.1860531159,-2.5030047259,21.419793629,-65.255885330,83.318680481,-33.746922930]
const a2 = [-0.0906148351,0.4527842806,0.5962700728,-1.7241829131,-4.1302112531,13.776631870,-8.6728470368]
const b0 = [0.7240946941,2.2382791861,-4.0025849485,-21.003576815,26.855641363,206.55133841,-355.60235612]
const b1 = [-0.5755498075,0.6995095521,3.8925673390,-17.215471648,192.67226447,-161.82646165,-165.20769346]
const b2 = [0.0976883116,-0.2557574982,-9.1558561530,20.642075974,-38.804430052,93.626774077,-29.666905585]

"""
    mutable struct Mix

gfhgfh
"""
mutable struct Mix
    n::Int
    componentes::Vector{String}
    MW::typeof([1.0]u"g/mol")
    x::typeof([1.0])
    y::typeof([1.0])
    m::typeof([1.0])
    σ::typeof([1.0]u"Å")
    ϵ::typeof([1.0]u"K")
    σᵢⱼ::typeof([1.0 1.0;1.0 1.0]u"Å")
    ϵᵢⱼ::typeof([1.0 1.0;1.0 1.0]u"K")
    k::Array{Float64,2}
end

function Mix(MW::typeof([1.0]u"g/mol"),x::typeof([1.0]),k::Array{Float64,2})
    n = length(MW)
    y = x
    componentes = ["Componente $i" for i=1:n]
    Mix(n,componentes,MW,x,y,k)
end

function Mix(componentes::Vector{String},MW::typeof([1.0]u"g/mol"),x::typeof([1.0]),k::Array{Float64,2})
    n = length(MW)
    y = x
    Mix(n,componentes,MW,x,y,k)
end

function Mix(MW::typeof([1.0]u"g/mol"),x::typeof([1.0]),y::typeof([1.0]),k::Array{Float64,2})
    n = length(MW)
    componentes = ["Componente $i" for i=1:n]
    Mix(n,componentes,MW,x,y,k)
end

function Mix(n::Int,MW::typeof([1.0]u"g/mol"),x::typeof([1.0]),k::Array{Float64,2})
    y = x
    componentes = ["Componente $i" for i=1:n]
    Mix(n,componentes,MW,x,y,k)
end

function Mix(n::Int,componentes::Vector{String},MW::typeof([1.0]u"g/mol"),x::typeof([1.0]),k::Array{Float64,2})
    y = x
    Mix(n,componentes,MW,x,y,k)
end

function Mix(componentes::Vector{String},MW::typeof([1.0]u"g/mol"),x::typeof([1.0]),y::typeof([1.0]),k::Array{Float64,2})
    n = length(MW) 
    Mix(n,componentes,MW,x,y,k)
end

function Mix(n::Int,MW::typeof([1.0]u"g/mol"),x::typeof([1.0]),y::typeof([1.0]),k::Array{Float64,2})
    componentes = ["Componente $i" for i=1:n] 
    Mix(n,componentes,MW,x,y,k)
end

function Mix(n::Int,componentes::Vector{String},MW::typeof([1.0]u"g/mol"),x::typeof([1.0]),y::typeof([1.0]),k::Array{Float64,2})
    σ(MW) = q01+(MW-MCH4)*q11/MW+(MW-MCH4)*(MW-2*MCH4)*q21/MW^2
    m(MW) =(q02 +(MW - MCH4)*q12/MW+(MW - MCH4)*(MW - 2*MCH4)*q22/MW^2)*MW
    ϵ(MW) = q03 +(MW - MCH4)*q13/MW+(MW - MCH4)*(MW - 2*MCH4)*q23/MW^2

    m = [m(MW[i]) for i=1:n]
    σ = [σ(MW[i]) for i=1:n]
    ϵ = [ϵ(MW[i]) for i=1:n]

    σᵢⱼ = Array{Float64}(undef,n,n)*u"Å"
    ϵᵢⱼ = Array{Float64}(undef,n,n)*u"K"
    for i=1:n
        for j=1:n
            σᵢⱼ[i,j] = 0.5*(σ[i]+σ[j])
            ϵᵢⱼ[i,j] = sqrt(ϵ[i]*ϵ[j])*(1-k[i,j])
        end
    end

    Mix(n,componentes,MW,x,y,m,σ,ϵ,σᵢⱼ,ϵᵢⱼ,k)
end

function Mix(MW::typeof([1.0]u"g/mol"),x::typeof([1.0]),m::typeof([1.0]),σ::typeof([1.0]u"Å"),ϵ::typeof([1.0]u"K"),k::Array{Float64,2})
    n = length(MW)
    y = x
    componentes = ["Componente $i" for i=1:n]
    Mix(n,componentes,MW,x,y,m,σ,ϵ,k)
end

function Mix(componentes::Vector{String},MW::typeof([1.0]u"g/mol"),x::typeof([1.0]),m::typeof([1.0]),σ::typeof([1.0]u"Å"),ϵ::typeof([1.0]u"K"),k::Array{Float64,2})
    n = length(MW)
    y = x
    Mix(n,componentes,MW,x,y,m,σ,ϵ,k)
end

function Mix(MW::typeof([1.0]u"g/mol"),x::typeof([1.0]),y::typeof([1.0]),m::typeof([1.0]),σ::typeof([1.0]u"Å"),ϵ::typeof([1.0]u"K"),k::Array{Float64,2})
    n = length(MW)
    componentes = ["Componente $i" for i=1:n]
    Mix(n,componentes,MW,x,y,m,σ,ϵ,k)
end

function Mix(n::Int,MW::typeof([1.0]u"g/mol"),x::typeof([1.0]),m::typeof([1.0]),σ::typeof([1.0]u"Å"),ϵ::typeof([1.0]u"K"),k::Array{Float64,2})
    y = x
    componentes = ["Componente $i" for i=1:n]
    Mix(n,componentes,MW,x,y,m,σ,ϵ,k)
end

function Mix(n::Int,componentes::Vector{String},MW::typeof([1.0]u"g/mol"),x::typeof([1.0]),m::typeof([1.0]),σ::typeof([1.0]u"Å"),ϵ::typeof([1.0]u"K"),k::Array{Float64,2})
    y = x
    Mix(n,componentes,MW,x,y,m,σ,ϵ,k)
end

function Mix(componentes::Vector{String},MW::typeof([1.0]u"g/mol"),x::typeof([1.0]),y::typeof([1.0]),m::typeof([1.0]),σ::typeof([1.0]u"Å"),ϵ::typeof([1.0]u"K"),k::Array{Float64,2})
    n = length(MW) 
    Mix(n,componentes,MW,x,y,m,σ,ϵ,k)
end

function Mix(n::Int,MW::typeof([1.0]u"g/mol"),x::typeof([1.0]),y::typeof([1.0]),m::typeof([1.0]),σ::typeof([1.0]u"Å"),ϵ::typeof([1.0]u"K"),k::Array{Float64,2})
    componentes = ["Componente $i" for i=1:n] 
    Mix(n,componentes,MW,x,y,m,σ,ϵ,k)
end

function Mix(n::Int,componentes::Vector{String},MW::typeof([1.0]u"g/mol"),x::typeof([1.0]),y::typeof([1.0]),m::typeof([1.0]),σ::typeof([1.0]u"Å"),ϵ::typeof([1.0]u"K"),k::Array{Float64,2})
    σᵢⱼ = Array{Float64}(undef,n,n)*u"Å"
    ϵᵢⱼ = Array{Float64}(undef,n,n)*u"K"
    for i=1:n
        for j=1:n
            σᵢⱼ[i,j] = 0.5*(σ[i]+σ[j])
            ϵᵢⱼ[i,j] = sqrt(ϵ[i]*ϵ[j])*(1-k[i,j])
        end
    end

    Mix(n,componentes,MW,x,y,m,σ,ϵ,σᵢⱼ,ϵᵢⱼ,k)
end

const mix = Mix

function zdisp(T,mix,ρ,phase)
    n = mix.n
    m,σ,ϵ,k = mix.m,mix.σ,mix.ϵ,mix.k
    σᵢⱼ,ϵᵢⱼ = mix.σᵢⱼ,mix.ϵᵢⱼ
    phase == "gas" ? x = mix.y : x = mix.x
    T = T |> u"K"

    d = σ.*(1 .-0.12*exp.(-3*ϵ/T)) 
    m̄ = ∑(x.*m)

    a = a0+ (m̄-1)*a1/m̄ + (m̄-1)*(m̄-2)*a2/m̄^2
    b = b0+ (m̄-1)*b1/m̄ + (m̄-1)*(m̄-2)*b2/m̄^2

    η = ∑(x.*m.*d.^3)
    η = η*π/6*ρ

    ∂ηI₁ = ∑(a[j]*j*η^(j-1) for j=1:7)
    ∂ηI₂ = ∑(b[j]*j*η^(j-1) for j=1:7)

    term1 = m̄*(8*η-2*η^2)/(1-η)^4;
    term2 = (1-m̄)*(20*η-27*η^2+12*η^3-2*η^4)/((1-η)*(2-η))^2;
    C₁ = (1+term1 + term2)^-1

    term1 = m̄*(-4*η^2+20*η+8)/(1-η)^5;
    term2 = (1-m̄)*(2*η^3+12*η^2-48*η+40)/((1-η)*(2-η))^3;
    C₂ = -C₁^2*(term1 + term2)

    m²ϵσ³ = 0u"Å^3"
    m²ϵ²σ³ = 0u"Å^3"

    for j=1:n
        for i=1:n
            m²ϵσ³ += x[i]*x[j]*m[i]*m[j]*(ϵᵢⱼ[i,j]/T)*σᵢⱼ[i,j]^3 
            m²ϵ²σ³ += x[i]*x[j]*m[i]*m[j]*(ϵᵢⱼ[i,j]/T)^2*σᵢⱼ[i,j]^3 
        end
    end

    I₂ = ∑(b[i]*η^(i-1) for i=1:7)

    z_disp = -2*π*ρ*∂ηI₁*m²ϵσ³ -π*ρ*m̄*(C₁*∂ηI₂+C₂*η*I₂)*m²ϵ²σ³

end

function z_hc(T,mix,ρ,phase)
    n = mix.n
    m,σ,ϵ,k = mix.m,mix.σ,mix.ϵ,mix.k
    phase == "gas" ? x = mix.y : x = mix.x
    T = T |> u"K"

    d = σ.*(1 .-0.12*exp.(-3*ϵ/T)) 
    m̄ = ∑(x.*m)

    ξ(nn) = (π*ρ/6)*∑(x[i]*m[i]*d[i]^nn for i=1:n)
    #ρ=ρ
    #η=η
    ξ = [ξ(i) for i=0:3]
    dᵢⱼ(i,j) = d[i]*d[j]/(d[i]+d[j])

    gʰˢ = Array{Float64}(undef,n,n)
    for j=1:n
        for i=1:n
            term1 = 1/(1-ξ[4]);
            term2 = dᵢⱼ(i,j)*3*ξ[3]/(1-ξ[4])^2;
            term3 = (dᵢⱼ(i,j))^2*2*ξ[3]^2/(1-ξ[4])^3;
            gʰˢ[i,j] = term1 + term2 + term3
        end
    end

    term1 = ξ[4]/(1-ξ[4]);
    term2 = 3*ξ[2]*ξ[3]/(ξ[1]*(1-ξ[4])^2);
    term3 = (3*ξ[3]^3-ξ[4]*ξ[3]^3)/(ξ[1]*(1-ξ[4])^3);
    Zʰˢ = term1 + term2 + term3

    ∂g = Array{Float64}(undef,n,n)
    for j=1:n
        for i=1:n
            term1 = ξ[4]/(1-ξ[4])^2;
            term2 = dᵢⱼ(i,j)*(3*ξ[3]/(1-ξ[4])^2+6*ξ[3]*ξ[4]/(1-ξ[4])^3);
            term3 = (dᵢⱼ(i,j))^2*(4*ξ[3]^2/(1-ξ[4])^3 + 6*ξ[3]^2*ξ[4]/(1-ξ[4])^4);
            ∂g[i,j] = term1 + term2 + term3
        end
    end

    z_hc = m̄*Zʰˢ - ∑(x[i]*(m[i]-1)*∂g[i,i]*gʰˢ[i,i]^-1 for i=1:n)

end

function obj_saft(η,T,P,mix,phase)
    n = mix.n
    m,σ,ϵ,k = mix.m,mix.σ,mix.ϵ,mix.k
    phase == "gas" ? x = mix.y : x = mix.x
    T = T |> u"K"
    P = P |> u"Pa"

    d = σ.*(1 .-0.12*exp.(-3*ϵ/T)) 
    m̄ = ∑(x.*m)

    η_new = ∑(x.*m.*d.^3)
    ρ = 6/π*η*η_new^-1
    Zʰᶜ = z_hc(T,mix,ρ,phase)
    Zᵈⁱˢᵖ = zdisp(T,mix,ρ,phase)

    Zcalc = 1 + Zʰᶜ + Zᵈⁱˢᵖ
    Pcalc = Zcalc*T*ρ*1.380649e-23u"J/K"
    Pcalc = Pcalc |> u"Pa"

    ustrip(abs(P-Pcalc))
end

function compresivilidad(η,T,mix,phase)
    n = mix.n
    m,σ,ϵ,k = mix.m,mix.σ,mix.ϵ,mix.k
    phase == "gas" ? x = mix.y : x = mix.x
    T = T |> u"K"

    d = σ.*(1 .-0.12*exp.(-3*ϵ/T)) 
    m̄ = ∑(x.*m)

    η_new = ∑(x.*m.*d.^3)
    ρ = 6/π*η*η_new^-1
    Zʰᶜ = z_hc(T,mix,ρ,phase)
    Zᵈⁱˢᵖ = zdisp(T,mix,ρ,phase)

    Zcalc = 1 + Zʰᶜ + Zᵈⁱˢᵖ
end

function compr(T,P,mix,phase)
    n = mix.n
    m,σ,ϵ,k = mix.m,mix.σ,mix.ϵ,mix.k
    x,M,y = mix.x,mix.MW,mix.y
    T = T |> u"K"

    if phase == "gas"
        ρ = P/(0.8*(R/Nₐ)*T) |> u"Å^-3"
        ηini = ∑(y.*m.*σ.^3)
        ηini = ηini*π/6*ρ
    elseif phase == "liq"
        MW = ∑(M)/n
        ρ = (800000u"g/m^3")*Nₐ/MW |> u"Å^-3"
        ηini = ∑(x.*m.*σ.^3)
        ηini = ηini*π/6*ρ
    else 
        error("No especifico el tipo de fase")
    end
    obj_saft1(η) = obj_saft(abs(η),T,P,mix,phase)
    ηsol = nlsolve(n_ary(obj_saft1),[ηini],iterations=100)
    η = ηsol.zero
    #println("P, $P, x, $x, η, $η, $phase")
    z = compresivilidad(η...,T,mix,phase)
    return z,η
end

function fugF(T,P,mix,phase)
    T = T |> u"K"
  
    Z,η = compr(T,P,mix,phase)
    ρ = P/(Z*kᵦ*T) |> u"Å^-3"
    μʰᶜ = mu_hc(T,mix,ρ,phase)
    μᵈⁱˢᵖ = mu_disp(T,mix,ρ,phase)
    μ = exp.(μʰᶜ + μᵈⁱˢᵖ .-log(Z))
    return μ
end

function mu_hc(T,mix,ρ,phase)
    n = mix.n
    m,σ,ϵ,k = mix.m,mix.σ,mix.ϵ,mix.k
    phase == "gas" ? x = mix.y : x = mix.x
    T = T |> u"K"

    d = σ.*(1 .-0.12*exp.(-3*ϵ/T)) 
    m̄ = ∑(x.*m)

    ξ(nn) = (π*ρ/6)*∑(x[i]*m[i]*d[i]^nn for i=1:n)
    #ρ=ρ
    #η=η
    ξ = [ξ(i) for i=0:3]
    dᵢⱼ(i,j) = d[i]*d[j]/(d[i]+d[j])
    dᵢⱼ = d.*d./(d+d)

    gʰˢ = Vector{Float64}(undef,n)
    for i=1:n
        term1 = 1/(1-ξ[4]);
        term2 = dᵢⱼ[i]*3*ξ[3]/(1-ξ[4])^2;
        term3 = (dᵢⱼ[i])^2*2*ξ[3]^2/(1-ξ[4])^3;
        gʰˢ[i] = term1 + term2 + term3
    end


    Zʰᶜ = z_hc(T,mix,ρ,phase)

    term1 = 3*ξ[2]*ξ[3]/(1-ξ[4]);
    term2 = ξ[3]^3/(ξ[4]*(1-ξ[4])^2);
    term3 = (ξ[3]^3/ξ[4]^2-ξ[1])*log(1-ξ[4]);
    aʰˢ = (1/ξ[1])*(term1 + term2 + term3)

    Aʰᶜ = m̄*aʰˢ-∑(x[i]*(m[i]-1)*log(gʰˢ[i]) for i=1:n)

    ξₖ = [(π*ρ/6)*m[i]*d[i]^(j-1) for j=1:4,i=1:n]

    ∂aʰˢ∂x = Vector{Float64}(undef,n)
    for i=1:n
        term1 = -ξₖ[1,i]/ξ[1]*aʰˢ;
        term2 = 3*(ξₖ[2,i]*ξ[3] + ξ[2]*ξₖ[3,i])/(1-ξ[4]);
        term3 = 3*ξ[2]*ξ[3]*ξₖ[4,i]/(1-ξ[4])^2;
        term4 = 3*ξ[3]^2*ξₖ[3,i]/(ξ[4]*(1-ξ[4])^2);
        term5 = ξ[3]^3*ξₖ[4,i]*(3*ξ[4]-1)/(ξ[4]^2*(1-ξ[4])^3);
        term6 = ((3*ξ[3]^2*ξₖ[3,i]*ξ[4]-2*ξ[3]^3*ξₖ[4,i])/ξ[4]^3-ξₖ[1,i])*log(1-ξ[4]);
        term7 = (ξ[1]-ξ[3]^3/ξ[4]^2)*ξₖ[4,i]/(1-ξ[4]);
        
        ∂aʰˢ∂x[i] = term1 + 1/ξ[1]*(term2 + term3 + term4 + term5 + term6 +term7)
    end

    ∂g∂x = Array{Float64}(undef,n,n)
    for k=1:n
        for i=1:n
            term1 = ξₖ[4,k]/(1-ξ[4])^2;
            term2 = dᵢⱼ[i]*(3*ξₖ[3,k]/(1-ξ[4])^2 + 6*ξ[3]*ξₖ[4,k]/(1-ξ[4])^3);
            term3 = (dᵢⱼ[i])^2*(4*ξ[3]*ξₖ[3,k]/(1-ξ[4])^3 + 6*ξ[3]^2*ξₖ[4,k]/(1-ξ[4])^4);
            ∂g∂x[i,k] = term1 + term2 + term3
        end
    end

    ∂aʰᶜ∂x = m*aʰˢ + m̄*∂aʰˢ∂x -(m.-1).*log.(gʰˢ) -∑(x[i]*(m[i]-1)*∂g∂x[i,:]*gʰˢ[i]^(-1) for i=1:n)

    μ = Aʰᶜ + Zʰᶜ .+ ∂aʰᶜ∂x .- ∑(x[i]*∂aʰᶜ∂x[i] for i=1:n)
end

function HelmholtzDisp(T,mix,ρ,phase)
    n = mix.n
    m,σ,ϵ,k = mix.m,mix.σ,mix.ϵ,mix.k
    σᵢⱼ,ϵᵢⱼ = mix.σᵢⱼ,mix.ϵᵢⱼ
    phase == "gas" ? x = mix.y : x = mix.x
    T = T |> u"K"

    d = σ.*(1 .-0.12*exp.(-3*ϵ/T)) 
    m̄ = ∑(x.*m)

    a = a0+ (m̄-1)*a1/m̄ + (m̄-1)*(m̄-2)*a2/m̄^2
    b = b0+ (m̄-1)*b1/m̄ + (m̄-1)*(m̄-2)*b2/m̄^2

    η = ∑(x.*m.*d.^3)
    η = η*π/6*ρ
    
    term1 = m̄*(8*η-2*η^2)/(1-η)^4;
    term2 = (1-m̄)*(20*η-27*η^2+12*η^3-2*η^4)/((1-η)*(2-η))^2;
    C₁ = (1+term1 + term2)^-1

    m²ϵσ³ = 0u"Å^3"
    m²ϵ²σ³ = 0u"Å^3"

    for j=1:n
        for i=1:n
            m²ϵσ³ += x[i]*x[j]*m[i]*m[j]*(ϵᵢⱼ[i,j]/T)*σᵢⱼ[i,j]^3 
            m²ϵ²σ³ += x[i]*x[j]*m[i]*m[j]*(ϵᵢⱼ[i,j]/T)^2*σᵢⱼ[i,j]^3 
        end
    end

    I₁ = ∑(a[i]*η^(i-1) for i=1:7)
    I₂ = ∑(b[i]*η^(i-1) for i=1:7)

    A_disp = -2*π*ρ*I₁*m²ϵσ³-π*ρ*m̄*C₁*I₂*m²ϵ²σ³

end

function mu_disp(T,mix,ρ,phase)
    n = mix.n
    m,σ,ϵ,k = mix.m,mix.σ,mix.ϵ,mix.k
    σᵢⱼ,ϵᵢⱼ = mix.σᵢⱼ,mix.ϵᵢⱼ
    phase == "gas" ? x = mix.y : x = mix.x
    T = T |> u"K"

    d = σ.*(1 .-0.12*exp.(-3*ϵ/T)) 
    m̄ = ∑(x.*m)

    a = a0+ (m̄-1)*a1/m̄ + (m̄-1)*(m̄-2)*a2/m̄^2
    b = b0+ (m̄-1)*b1/m̄ + (m̄-1)*(m̄-2)*b2/m̄^2

    η = ∑(x.*m.*d.^3)
    η = η*π/6*ρ

    Zᵈⁱˢᵖ = zdisp(T,mix,ρ,phase)
    Aᵈⁱˢᵖ = HelmholtzDisp(T,mix,ρ,phase)

    ξₖ = [(π*ρ/6)*m[i]*d[i]^(j-1) for j=1:4,i=1:n]

    m²ϵσ³ = 0u"Å^3"
    m²ϵ²σ³ = 0u"Å^3"

    for j=1:n
        for i=1:n
            m²ϵσ³ += x[i]*x[j]*m[i]*m[j]*(ϵᵢⱼ[i,j]/T)*σᵢⱼ[i,j]^3 
            m²ϵ²σ³ += x[i]*x[j]*m[i]*m[j]*(ϵᵢⱼ[i,j]/T)^2*σᵢⱼ[i,j]^3 
        end
    end

    m²ϵσ³ₓₖ = 2*m.*∑(x[j]*m[j]*(ϵᵢⱼ[:,j]/T).*σᵢⱼ[:,j].^3 for j=1:n)
    m²ϵ²σ³ₓₖ = 2*m.*∑(x[j]*m[j]*((ϵᵢⱼ[:,j]/T).^2).*σᵢⱼ[:,j].^3 for j=1:n)

    I₁ = ∑(a[i]*η^(i-1) for i=1:7)
    I₂ = ∑(b[i]*η^(i-1) for i=1:7)

    term1 = m̄*(8*η-2*η^2)/(1-η)^4;
    term2 = (1-m̄)*(20*η-27*η^2+12*η^3-2*η^4)/((1-η)*(2-η))^2;
    C₁ = (1+term1 + term2)^-1

    term1 = m̄*(-4*η^2+20*η+8)/(1-η)^5;
    term2 = (1-m̄)*(2*η^3+12*η^2-48*η+40)/((1-η)*(2-η))^3;
    C₂ = -C₁^2*(term1 + term2)

    C₁ₓₖ = C₂*ξₖ[4,:] - (C₁^2)*(m*(8*η-2*η^2)/(1-η)^4 -m*(20*η-27*η^2+12*η^3-2*η^4)/((1-η)*(2-η))^2)
    
    aₖ = a1*m'/m̄^2 +a2*m'*(3-4/m̄)/m̄^2
    bₖ = b1*m'/m̄^2 +b2*m'*(3-4/m̄)/m̄^2

    I₁ₓₖ = Vector{Float64}(undef,n)
    I₂ₓₖ = Vector{Float64}(undef,n)
    for j=1:n
        sum1 = ∑((a[i]*(i-1)*ξₖ[4,j]*η^(i-2))+aₖ[i,j]*η^(i-1) for i=1:7)
        sum2 = ∑((b[i]*(i-1)*ξₖ[4,j]*η^(i-2))+bₖ[i,j]*η^(i-1) for i=1:7)
        I₁ₓₖ[j] = sum1
        I₂ₓₖ[j] = sum2
    end

    ∂aᵈⁱˢᵖ∂x = -2*π*ρ*(I₁ₓₖ*m²ϵσ³+I₁*m²ϵσ³ₓₖ) -π*ρ*((m*C₁*I₂+m̄*C₁ₓₖ*I₂+m̄*C₁*I₂ₓₖ)*m²ϵ²σ³ +m̄*C₁*I₂*m²ϵ²σ³ₓₖ)

    μ = Aᵈⁱˢᵖ + Zᵈⁱˢᵖ .+ ∂aᵈⁱˢᵖ∂x .- ∑(x[i]*∂aᵈⁱˢᵖ∂x[i] for i=1:n)
    
end

function obj_buble(T,P,mix)
    iter = 0
    error = 100
    K = Vector{Float64}(undef, mix.n)
    while error > 1e-5 && iter < 100
        iter +=1
        fliq = fugF(T,P,mix,"liq")
        fgas = fugF(T,P,mix,"gas")
        K = fliq./fgas
        y_calc = mix.x.*K
        y_sum = ∑(y_calc)
        #println("mix.y= $(mix.y), y_calc= $y_calc")
        error = ∑(abs.(y_calc-mix.y))
        mix.y = y_calc/y_sum
        P = y_sum*P
    end

    return P,mix.y

end

"""
    pxy(T,P,mix::mix,xx;uni=u"Torr")

**Los campos de entrada son:**
- `T :: Float` Es la temperatura del sistema
- `P :: m-element Array{Float64,1}` presión inicial
- `mix` Datos de entrada de la mezcla
- `xx :: StepRange ó Vector` Son los puntos en el liquido donde se buscara el equilibrio con el vapor
- `uni` Son las unidades en las que se desea el resultado, por defecto esta en Torr

**Salida:**

Regresa las fracciones de liquido `x`,
las fracciones del vapor `y` y la presión `P` en Torr (`m` es la longitud de `xx`).

- `x :: n×m Array{Float64,2}` liquido
- `y :: n×m Array{Float64,2}` vapor
- `P :: m-element Array{Float64,1}` presión
"""
function pxy(T,P,mix::mix,xx;uni=u"Torr")
    T = T |> u"K"
    P = P |> u"Torr"
    xx = general.nc(xx...)
    m = size(xx,2)
    y = Array{Float64}(undef,mix.n,m)
    p = Vector{Float64}(undef,m)*u"Torr"

    mix.x = xx[:,1]
    mix.y = mix.x 
    p[1],y[:,1] = obj_buble(T,P,mix)
    for i=2:m
        mix.x = xx[:,i]
        mix.y = y[:,i-1] 
        p[i],y[:,i] = obj_buble(T,p[i-1],mix)
    end

    p = p .|> uparse(uni)
    return xx,y,p
    
end

end