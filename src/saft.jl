using Unitful
using Plots
using NLsolve

Unitful.preferunits(u"K",u"g",u"mol",u"Å")

upreferred([1u"kg",1u"g",1u"cg",1u"mg"])

∑(x)=sum(x)

1+11

mutable struct Mix
    n::Int
    componentes::Vector{String}
    MW::typeof([1.0]u"g/mol")
    x::typeof([1.0])
    y::typeof([1.0])
    m::typeof([1.0])
    σ::typeof([1.0]u"Å")
    ϵ::typeof([1.0]u"K")
    k::Array{Float64,2}
end

const mix = Mix

const q01, q11, q21 = [3.7039,-0.3226,0.6907]u"Å"
const q02, q12, q22 = [0.06233,-0.02236,-0.01563]u"mol/g"
const q03, q13, q23 = [150.03,80.68,38.96]u"K"
const MCH4 = 16.043u"g/mol"
const R = 8.31446261815324u"MPa*cm^3/(K*mol)"
const Nₐ = 6.02214076e23u"1/(mol)"
const kᵦ = R/Nₐ

function hidrocarburos(mix)
    MW = mix.MW

    σ(MW) = q01+(MW-MCH4)*q11/MW+(MW-MCH4)*(MW-2*MCH4)*q21/MW^2
    m(MW) =(q02 +(MW - MCH4)*q12/MW+(MW - MCH4)*(MW - 2*MCH4)*q22/MW^2)*MW
    ϵ(MW) = q03 +(MW - MCH4)*q13/MW+(MW - MCH4)*(MW - 2*MCH4)*q23/MW^2

    m = [m(MW[i]) for i=1:mix.n]
    σ = [σ(MW[i]) for i=1:mix.n]
    ϵ = [ϵ(MW[i]) for i=1:mix.n]

    return m,σ,ϵ

end

mix(2,["metano","butano"],[16.043,58.123]u"kg/mol",[5,.5],[.1,.9],[.1,.11],[.1,.9]u"Å",[.1,.9]u"K",[0.0 0.022;0.022 0.0])
mix1 = mix(2,["metano","butano"],[16.043,58.123]u"kg/mol",[.5,.5],[.1,.9],[.1,.11],[.1,.9]u"Å",[.1,.9]u"K",[0.0 0.022;0.022 0.0])
mix1.MW
mix1.MW = [16.043,58.123]u"g/mol"

mix1.m,mix1.σ,mix1.ϵ = hidrocarburos(mix1)

const a0 = [0.9105631445,0.6361281449,2.6861347891,-26.547362491,97.759208784,-159.59154087,91.297774084]
const a1 = [-0.3084016918,0.1860531159,-2.5030047259,21.419793629,-65.255885330,83.318680481,-33.746922930]
const a2 = [-0.0906148351,0.4527842806,0.5962700728,-1.7241829131,-4.1302112531,13.776631870,-8.6728470368]
const b0 = [0.7240946941,2.2382791861,-4.0025849485,-21.003576815,26.855641363,206.55133841,-355.60235612]
const b1 = [-0.5755498075,0.6995095521,3.8925673390,-17.215471648,192.67226447,-161.82646165,-165.20769346]
const b2 = [0.0976883116,-0.2557574982,-9.1558561530,20.642075974,-38.804430052,93.626774077,-29.666905585]

function zdisp(T,M,ρ,x,k)
    n = length(x)
    T = T |> u"K"

    σ(M) = q01+(M-MCH4)*q11/M+(M-MCH4)*(M-2*MCH4)*q21/M^2
    m(M) =(q02 +(M - MCH4)*q12/M+(M - MCH4)*(M - 2*MCH4)*q22/M^2)*M
    ϵ(M) = q03 +(M - MCH4)*q13/M+(M - MCH4)*(M - 2*MCH4)*q23/M^2

    m = [m(M[i]) for i=1:n]
    σ = [σ(M[i]) for i=1:n]
    ϵ = [ϵ(M[i]) for i=1:n]

    d = σ.*(1 .-0.12*exp.(-3*ϵ/T)) 
    m̄ = ∑(x.*m)

    a = a0+ (m̄-1)*a1/m̄ + (m̄-1)*(m̄-2)*a2/m̄^2
    b = b0+ (m̄-1)*b1/m̄ + (m̄-1)*(m̄-2)*b2/m̄^2

    η = ∑(x.*m.*d.^3)
    η = η*π/6*ρ

    σᵢⱼ = Array{Float64}(undef,n,n)*u"Å"
    ϵᵢⱼ = Array{Float64}(undef,n,n)*u"K"
    for i=1:n
        for j=1:n
            σᵢⱼ[i,j] = 0.5*(σ[i]+σ[j])
            ϵᵢⱼ[i,j] = sqrt(ϵ[i]*ϵ[j])*(1-k[i,j])
        end
    end

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

    for i=1:n
        for j=1:n
            m²ϵσ³ += x[i]*x[j]*m[i]*m[j]*(ϵᵢⱼ[i,j]/T)*σᵢⱼ[i,j]^3 
            m²ϵ²σ³ += x[i]*x[j]*m[i]*m[j]*(ϵᵢⱼ[i,j]/T)^2*σᵢⱼ[i,j]^3 
        end
    end

    I₂ = ∑(b[i]*η^(i-1) for i=1:7)

    z_disp = -2*π*ρ*∂ηI₁*m²ϵσ³ -π*ρ*m̄*(C₁*∂ηI₂+C₂*η*I₂)*m²ϵ²σ³

end

function z_hc(T,M,ρ,x)
    n = length(x)
    T = T |> u"K"

    σ(M) = q01+(M-MCH4)*q11/M+(M-MCH4)*(M-2*MCH4)*q21/M^2
    m(M) =(q02 +(M - MCH4)*q12/M+(M - MCH4)*(M - 2*MCH4)*q22/M^2)*M
    ϵ(M) = q03 +(M - MCH4)*q13/M+(M - MCH4)*(M - 2*MCH4)*q23/M^2

    m = [m(M[i]) for i=1:n]
    σ = [σ(M[i]) for i=1:n]
    ϵ = [ϵ(M[i]) for i=1:n]

    d = σ.*(1 .-0.12*exp.(-3*ϵ/T)) 
    m̄ = ∑(x.*m)

    ξ(nn) = (π*ρ/6)*∑(x[i]*m[i]*d[i]^nn for i=1:n)
    #ρ=ρ
    #η=η
    ξ = [ξ(i) for i=0:3]
    dᵢⱼ(i,j) = d[i]*d[j]/(d[i]+d[j])

    gʰˢ = Array{Float64}(undef,n,n)
    for i=1:n
        for j=1:n
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
    for i=1:n
        for j=1:n
            term1 = ξ[4]/(1-ξ[4])^2;
            term2 = dᵢⱼ(i,j)*(3*ξ[3]/(1-ξ[4])^2+6*ξ[3]*ξ[4]/(1-ξ[4])^3);
            term3 = (dᵢⱼ(i,j))^2*(4*ξ[3]^2/(1-ξ[4])^3 + 6*ξ[3]^2*ξ[4]/(1-ξ[4])^4);
            ∂g[i,j] = term1 + term2 + term3
        end
    end

    z_hc = m̄*Zʰˢ - ∑(x[i]*(m[i]-1)*∂g[i,i]*gʰˢ[i,i]^-1 for i=1:n)

end

function obj_saft(η,T,P,M,x,k)
    n = length(x)
    T = T |> u"K"
    P = P |> u"Pa"

    σ(M) = q01+(M-MCH4)*q11/M+(M-MCH4)*(M-2*MCH4)*q21/M^2
    m(M) =(q02 +(M - MCH4)*q12/M+(M - MCH4)*(M - 2*MCH4)*q22/M^2)*M
    ϵ(M) = q03 +(M - MCH4)*q13/M+(M - MCH4)*(M - 2*MCH4)*q23/M^2

    m = [m(M[i]) for i=1:n]
    σ = [σ(M[i]) for i=1:n]
    ϵ = [ϵ(M[i]) for i=1:n]

    d = σ.*(1 .-0.12*exp.(-3*ϵ/T)) 
    m̄ = ∑(x.*m)

    η_new = ∑(x.*m.*d.^3)
    ρ = 6/π*η*η_new^-1
    Zʰᶜ = z_hc(T,M,ρ,x)
    Zᵈⁱˢᵖ = zdisp(T,M,ρ,x,k)

    Zcalc = 1 + Zʰᶜ + Zᵈⁱˢᵖ
    Pcalc = Zcalc*T*ρ*1.380649e-23u"J/K"
    Pcalc = Pcalc |> u"Pa"

    ustrip(abs(P-Pcalc))
end

function compresivilidad(η,T,M,x,k)
    n = length(x)
    T = T |> u"K"

    σ(M) = q01+(M-MCH4)*q11/M+(M-MCH4)*(M-2*MCH4)*q21/M^2
    m(M) =(q02 +(M - MCH4)*q12/M+(M - MCH4)*(M - 2*MCH4)*q22/M^2)*M
    ϵ(M) = q03 +(M - MCH4)*q13/M+(M - MCH4)*(M - 2*MCH4)*q23/M^2

    m = [m(M[i]) for i=1:n]
    σ = [σ(M[i]) for i=1:n]
    ϵ = [ϵ(M[i]) for i=1:n]

    d = σ.*(1 .-0.12*exp.(-3*ϵ/T)) 
    m̄ = ∑(x.*m)

    η_new = ∑(x.*m.*d.^3)
    ρ = 6/π*η*η_new^-1
    Zʰᶜ = z_hc(T,M,ρ,x)
    Zᵈⁱˢᵖ = zdisp(T,M,ρ,x,k)

    Zcalc = 1 + Zʰᶜ + Zᵈⁱˢᵖ
end

function compr(T,P,M,x,k,phase)
    n = length(x)
    T = T |> u"K"
    σ(M) = q01+(M-MCH4)*q11/M+(M-MCH4)*(M-2*MCH4)*q21/M^2
    m(M) =(q02 +(M - MCH4)*q12/M+(M - MCH4)*(M - 2*MCH4)*q22/M^2)*M

    m = [m(M[i]) for i=1:n]
    σ = [σ(M[i]) for i=1:n]

    if phase == "gas"
        ρ = P/(0.8*(R/Nₐ)*T) |> u"Å^-3"
        ηini = ∑(x.*m.*σ.^3)
        ηini = ηini*π/6*ρ
    elseif phase == "liq"
        MW = ∑(M)/n
        ρ = (800000u"g/m^3")*Nₐ/MW |> u"Å^-3"
        ηini = ∑(x.*m.*σ.^3)
        ηini = ηini*π/6*ρ
    else 
        error("No especifico el tipo de fase")
    end
    obj_saft1(η) = obj_saft(abs(η),T,P,M,x,k)
    ηsol = nlsolve(n_ary(obj_saft1),[ηini])
    η = ηsol.zero
    println("P, $P, x, $x, η, $η, $phase")
    z = compresivilidad(η...,T,M,x,k)
    return z,η
end

function fugF(T,P,M,x,k,phase)
    n = length(x)
    T = T |> u"K"
    Z,η = compr(T,P,M,x,k,phase)
    ρ = P/(Z*kᵦ*T) |> u"Å^-3"
    μʰᶜ = mu_hc(T,M,ρ,x)
    μᵈⁱˢᵖ = mu_disp(T,M,ρ,x,k)
    μ = exp.(μʰᶜ + μᵈⁱˢᵖ .-log(Z))
    #return μʰᶜ,μᵈⁱˢᵖ
end

function mu_hc(T,M,ρ,x)
    n = length(x)
    T = T |> u"K"

    σ(M) = q01+(M-MCH4)*q11/M+(M-MCH4)*(M-2*MCH4)*q21/M^2
    m(M) =(q02 +(M - MCH4)*q12/M+(M - MCH4)*(M - 2*MCH4)*q22/M^2)*M
    ϵ(M) = q03 +(M - MCH4)*q13/M+(M - MCH4)*(M - 2*MCH4)*q23/M^2

    m = [m(M[i]) for i=1:n]
    σ = [σ(M[i]) for i=1:n]
    ϵ = [ϵ(M[i]) for i=1:n]

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


    Zʰᶜ = z_hc(T,M,ρ,x)

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
    for i=1:n
        for k=1:n
            term1 = ξₖ[4,k]/(1-ξ[4])^2;
            term2 = dᵢⱼ[i]*(3*ξₖ[3,k]/(1-ξ[4])^2 + 6*ξ[3]*ξₖ[4,k]/(1-ξ[4])^3);
            term3 = (dᵢⱼ[i])^2*(4*ξ[3]*ξₖ[3,k]/(1-ξ[4])^3 + 6*ξ[3]^2*ξₖ[4,k]/(1-ξ[4])^4);
            ∂g∂x[i,k] = term1 + term2 + term3
        end
    end

    ∂aʰᶜ∂x = m*aʰˢ + m̄*∂aʰˢ∂x -(m.-1).*log.(gʰˢ) -∑(x[i]*(m[i]-1)*∂g∂x[i,:]*gʰˢ[i]^(-1) for i=1:n)

    μ = Aʰᶜ + Zʰᶜ .+ ∂aʰᶜ∂x .- ∑(x[i]*∂aʰᶜ∂x[i] for i=1:n)

end

function HelmholtzDisp(T,M,ρ,x,k)
    n = length(x)
    T = T |> u"K"
    a0 = Vector{Float64}(undef,7)
    a1 = Vector{Float64}(undef,7)
    a2 = Vector{Float64}(undef,7)
    b0 = Vector{Float64}(undef,7)
    b1 = Vector{Float64}(undef,7)
    b2 = Vector{Float64}(undef,7)
    a0[1] = 0.9105631445;
    a0[2] = 0.6361281449;
    a0[3] = 2.6861347891;
    a0[4] = -26.547362491;
    a0[5] = 97.759208784;
    a0[6] = -159.59154087;
    a0[7] = 91.297774084;
    a1[1] = -0.3084016918;
    a1[2] = 0.1860531159;
    a1[3] = -2.5030047259;
    a1[4] = 21.419793629;
    a1[5] = -65.255885330;
    a1[6] = 83.318680481;
    a1[7] = -33.746922930;
    a2[1] = -0.0906148351;
    a2[2] = 0.4527842806;
    a2[3] = 0.5962700728;
    a2[4] = -1.7241829131;
    a2[5] = -4.1302112531;
    a2[6] = 13.776631870;
    a2[7] = -8.6728470368;
    b0[1] = 0.7240946941;
    b0[2] = 2.2382791861;
    b0[3] = -4.0025849485;
    b0[4] = -21.003576815;
    b0[5] = 26.855641363;
    b0[6] = 206.55133841;
    b0[7] = -355.60235612;
    b1[1] = -0.5755498075;
    b1[2] = 0.6995095521;
    b1[3] = 3.8925673390;
    b1[4] = -17.215471648;
    b1[5] = 192.67226447;
    b1[6] = -161.82646165;
    b1[7] = -165.20769346;
    b2[1] = 0.0976883116;
    b2[2] = -0.2557574982;
    b2[3] = -9.1558561530;
    b2[4] = 20.642075974;
    b2[5] = -38.804430052;
    b2[6] = 93.626774077;
    b2[7] = -29.666905585;

    σ(M) = q01+(M-MCH4)*q11/M+(M-MCH4)*(M-2*MCH4)*q21/M^2
    m(M) =(q02 +(M - MCH4)*q12/M+(M - MCH4)*(M - 2*MCH4)*q22/M^2)*M
    ϵ(M) = q03 +(M - MCH4)*q13/M+(M - MCH4)*(M - 2*MCH4)*q23/M^2

    m = [m(M[i]) for i=1:n]
    σ = [σ(M[i]) for i=1:n]
    ϵ = [ϵ(M[i]) for i=1:n]

    d = σ.*(1 .-0.12*exp.(-3*ϵ/T)) 
    m̄ = ∑(x.*m)

    a = a0+ (m̄-1)*a1/m̄ + (m̄-1)*(m̄-2)*a2/m̄^2
    b = b0+ (m̄-1)*b1/m̄ + (m̄-1)*(m̄-2)*b2/m̄^2

    η = ∑(x.*m.*d.^3)
    η = η*π/6*ρ

    σᵢⱼ = Array{Float64}(undef,n,n)*u"Å"
    ϵᵢⱼ = Array{Float64}(undef,n,n)*u"K"
    for i=1:n
        for j=1:n
            σᵢⱼ[i,j] = 0.5*(σ[i]+σ[j])
            ϵᵢⱼ[i,j] = sqrt(ϵ[i]*ϵ[j])*(1-k[i,j])
        end
    end

    term1 = m̄*(8*η-2*η^2)/(1-η)^4;
    term2 = (1-m̄)*(20*η-27*η^2+12*η^3-2*η^4)/((1-η)*(2-η))^2;
    C₁ = (1+term1 + term2)^-1

    m²ϵσ³ = 0u"Å^3"
    m²ϵ²σ³ = 0u"Å^3"

    for i=1:n
        for j=1:n
            m²ϵσ³ += x[i]*x[j]*m[i]*m[j]*(ϵᵢⱼ[i,j]/T)*σᵢⱼ[i,j]^3 
            m²ϵ²σ³ += x[i]*x[j]*m[i]*m[j]*(ϵᵢⱼ[i,j]/T)^2*σᵢⱼ[i,j]^3 
        end
    end

    I₁ = ∑(a[i]*η^(i-1) for i=1:7)
    I₂ = ∑(b[i]*η^(i-1) for i=1:7)

    A_disp = -2*π*ρ*I₁*m²ϵσ³-π*ρ*m̄*C₁*I₂*m²ϵ²σ³

end

function mu_disp(T,M,ρ,x,k)
    n = length(x)
    T = T |> u"K"
    a0 = Vector{Float64}(undef,7)
    a1 = Vector{Float64}(undef,7)
    a2 = Vector{Float64}(undef,7)
    b0 = Vector{Float64}(undef,7)
    b1 = Vector{Float64}(undef,7)
    b2 = Vector{Float64}(undef,7)
    a0[1] = 0.9105631445;
    a0[2] = 0.6361281449;
    a0[3] = 2.6861347891;
    a0[4] = -26.547362491;
    a0[5] = 97.759208784;
    a0[6] = -159.59154087;
    a0[7] = 91.297774084;
    a1[1] = -0.3084016918;
    a1[2] = 0.1860531159;
    a1[3] = -2.5030047259;
    a1[4] = 21.419793629;
    a1[5] = -65.255885330;
    a1[6] = 83.318680481;
    a1[7] = -33.746922930;
    a2[1] = -0.0906148351;
    a2[2] = 0.4527842806;
    a2[3] = 0.5962700728;
    a2[4] = -1.7241829131;
    a2[5] = -4.1302112531;
    a2[6] = 13.776631870;
    a2[7] = -8.6728470368;
    b0[1] = 0.7240946941;
    b0[2] = 2.2382791861;
    b0[3] = -4.0025849485;
    b0[4] = -21.003576815;
    b0[5] = 26.855641363;
    b0[6] = 206.55133841;
    b0[7] = -355.60235612;
    b1[1] = -0.5755498075;
    b1[2] = 0.6995095521;
    b1[3] = 3.8925673390;
    b1[4] = -17.215471648;
    b1[5] = 192.67226447;
    b1[6] = -161.82646165;
    b1[7] = -165.20769346;
    b2[1] = 0.0976883116;
    b2[2] = -0.2557574982;
    b2[3] = -9.1558561530;
    b2[4] = 20.642075974;
    b2[5] = -38.804430052;
    b2[6] = 93.626774077;
    b2[7] = -29.666905585;

    σ(M) = q01+(M-MCH4)*q11/M+(M-MCH4)*(M-2*MCH4)*q21/M^2
    m(M) =(q02 +(M - MCH4)*q12/M+(M - MCH4)*(M - 2*MCH4)*q22/M^2)*M
    ϵ(M) = q03 +(M - MCH4)*q13/M+(M - MCH4)*(M - 2*MCH4)*q23/M^2

    m = [m(M[i]) for i=1:n]
    σ = [σ(M[i]) for i=1:n]
    ϵ = [ϵ(M[i]) for i=1:n]

    d = σ.*(1 .-0.12*exp.(-3*ϵ/T)) 
    m̄ = ∑(x.*m)

    a = a0+ (m̄-1)*a1/m̄ + (m̄-1)*(m̄-2)*a2/m̄^2
    b = b0+ (m̄-1)*b1/m̄ + (m̄-1)*(m̄-2)*b2/m̄^2

    η = ∑(x.*m.*d.^3)
    η = η*π/6*ρ

    σᵢⱼ = Array{Float64}(undef,n,n)*u"Å"
    ϵᵢⱼ = Array{Float64}(undef,n,n)*u"K"
    for i=1:n
        for j=1:n
            σᵢⱼ[i,j] = 0.5*(σ[i]+σ[j])
            ϵᵢⱼ[i,j] = sqrt(ϵ[i]*ϵ[j])*(1-k[i,j])
        end
    end

    Zᵈⁱˢᵖ = zdisp(T,M,ρ,x,k)
    Aᵈⁱˢᵖ = HelmholtzDisp(T,M,ρ,x,k)

    ξₖ = [(π*ρ/6)*m[i]*d[i]^(j-1) for j=1:4,i=1:n]

    m²ϵσ³ = 0u"Å^3"
    m²ϵ²σ³ = 0u"Å^3"

    for i=1:n
        for j=1:n
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

compr(21u"°C",100u"bar",[16.043,58.123]u"g/mol",[.5,.5],k,"gas")
compr(21u"°C",100u"bar",[16.043,58.123]u"g/mol",[.5,.5],k,"liq")
fugF(21u"°C",100u"bar",[16.043,58.123]u"g/mol",[.5,.5],k,"gas")
fugF(21u"°C",100u"bar",[16.043,58.123]u"g/mol",[.5,.5],k,"liq")
k=[0 0.022;
0.022 0]

zdisp(21u"°C",[16.043,58.123]u"g/mol",.01u"Å^-3",[.5,.5],k)
z_hc(21u"°C",[16.043,58.123]u"g/mol",.01u"Å^-3",[.5,.5])
obj_saft(0.097,21u"°C",100u"bar",[16.043,58.123]u"g/mol",[0.9085757578461996, 0.09142424215380038],k)
mu_hc(21u"°C",[16.043,58.123]u"g/mol",.01u"Å^-3",[.5,.5])
HelmholtzDisp(21u"°C",[16.043,58.123]u"g/mol",.01u"Å^-3",[.5,.5],k)
mu_disp(21u"°C",[16.043,58.123]u"g/mol",.01u"Å^-3",[.5,.5],k)

obj_saft1(η) = abs(obj_saft(abs(η),21u"°C",100u"bar",[16.043,58.123]u"g/mol",[.5,.5],k))
ηsol = nlsolve(n_ary(obj_saft1),[0.15],iterations=50000)
ηsol = nlsolve(n_ary(obj_saft1),[.7],method = :newton)
ηsol = nlsolve(n_ary(obj_saft1),[.5],method = :anderson)
ηsol = mcpsolve(n_ary(obj_saft1),[.1],[1.],[.55], reformulation = :smooth)
obj_saft(0.3024911711129972,21u"°C",100u"bar",[16.043,58.123]u"g/mol",[.5,.5],k)
η1 = ηsol.zero
converged(ηsol)


function f()
    foca = 5
    caca = 0
    ssds =89
end

[caca] = f()

r = mcpsolve(n_ary(f), [0.], [50.],[10.0],reformulation = :smooth, autodiff = :forward)

plot(xx->obj_saft(xx,21u"°C",100u"bar",[16.043,58.123]u"g/mol",[.5,.5],k),0.301,0.303)

function obj_buble(T,P,M,x,y,k)
    fliq = fugF(T,P,M,x,k,"liq")

    #yini = abs.(yini)
    #push!(yini,1-sum(yini))
    #ysum = sum(yini)
    #y = yini/ysum
    fgas = fugF(T,P,M,y,k,"gas")
    #display(y)
    error = sum(abs.(x.*fliq -y.*fgas))
end

function obj_y(T,P,M,fliq,x,y,k)
    val = 100.1
    iter = 0

    while val > 1e-5 && iter < 100
        iter =+ 1
        fgas = fugF(T,P,M,y,k,"gas")
        K = fliq./fgas
        y_calc = x.*K

        val = ∑(abs(y_old - y_calc))
        y_old = y_calc
        y = y_calc/∑(y_calc)
    end
    
    return y
end

function obj_buble(T,P,M,x,y,k)
    fliq = fugF(T,P,M,x,k,"liq")

    #yini = abs.(yini)
    #push!(yini,1-sum(yini))
    #ysum = sum(yini)
    #y = yini/ysum
    y = obj_y(T,P,M,fliq,x,y,k)
    display(y)
end

function obj_y(T,P,M,fliq,x,y,k)
    val = 100.1
    iter = 0
    y_old = y

    while val > 1e-5 && iter < 100
        iter =+ 1
        fgas = fugF(T,P,M,y,k,"gas")
        K = fliq./fgas
        y_calc = x.*K

        val = ∑(abs.(y_old - y_calc))
        y_old = y_calc
        println("y, $y")
        P = ∑(y_calc)*P
        y = y_calc/∑(y_calc)
        
    end
    
    return y
end

obj_buble1(y) = obj_buble(21u"°C",50u"bar",[16.043,58.123]u"g/mol",[.3,.7],y,k)
ηsol = nlsolve(obj_buble1,[.8;.2],iterations=50000)
obj_buble(21u"°C",50u"bar",[16.043,58.123]u"g/mol",[.3,.7],[.8,.2],k)

obj_buble1(p) = obj_buble(21u"°C",abs.(p)*u"bar",[16.043,58.123]u"g/mol",[.5,.5],[0.6,0.4],k)
ηsol = nlsolve(n_ary(obj_buble1),[125.],iterations=50)
obj_buble(21u"°C",100u"bar",[16.043,58.123]u"g/mol",[.5,.5],[0.6,0.4],k)
1+1
k=[0 3.00e-4 1.15e-2;
3.00e-4 0 5.10e-3;
1.15e-2 5.10e-3 0]

obj_buble1(p) = obj_buble(233.15u"K",abs.(p)*u"Pa",[16.043,30.07,44.096]u"g/mol",[.1,.3,.6],[0.6851,0.2442,0.0707],k)
ηsol = nlsolve(n_ary(obj_buble1),[1.6e6],iterations=50)
obj_buble(233.15u"K",1.0e6u"Pa",[16.043,30.07,44.096]u"g/mol",[.1,.3,.6],[0.727215096714321,0.1565049386734792,0.00038759300172423543],k)
obj_buble(233.15u"K",1.3e6u"Pa",[16.043,30.07,44.096]u"g/mol",[.1,.3,.6],[0.6851,0.2442,0.0707],k)
obj_buble(233.15u"K",1.3e6u"Pa",[16.043,30.07,44.096]u"g/mol",[.1,.3,.6],[2562.2888519723906, 891.721315066141, 231.31956888180963],k)


plot(xx->obj_buble(233.15u"K",xx*u"Pa",[16.043,30.07,44.096]u"g/mol",[.1,.3,.6],[0.6851,0.2442,0.0707],k),1e6,2e6)

[obj_buble(233.15u"K",xx*u"Pa",[16.043,30.07,44.096]u"g/mol",[.1,.3,.6],[0.6851,0.2442,0.0707],k) for xx=1e6:1e5:2e6]