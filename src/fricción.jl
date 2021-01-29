module fricción
    
export factor_friccion,laminar,Blasius,Moody,Altshul,Wood,Churchill,Colebrook

using NLsolve

laminar(Re) = 64/Re
Blasius(Re) = 0.3164/Re^0.25
Moody(Re,eD) = 1.375e-3*(1+(2e4*eD+1e6/Re)^(1/3))
Altshul(Re,eD) = 0.11*(68/Re +eD)^0.25
Wood(Re,eD) = 0.094*eD^0.225 +0.53*eD+88*(eD^.4)*Re^-(1.62*eD^.134)
function Churchill(Re,eD)
    fsol = nlsolve(n_ary((f) ->-2*log10(eD/3.7 +(7/Re)^.9) -1/√f),[.01])
    f = fsol.zero
    f[1]
end
function Colebrook(Re,eD)
    fsol = nlsolve(n_ary((f) ->-2*log10(eD/3.7 +2.51/(Re*√f)) -1/√f),[.01])
    f = fsol.zero
    f[1]
end

métodos = Dict(     "laminar"=> (0, 2300, nothing, nothing),
                    "Blasius"=> (3000, 200000, nothing, nothing),
                    "Colebrook"=> (-Inf, Inf, -Inf, Inf),
                    "Moody"=> (4e3, 1e8, 0, 1e-2),
                    "Altshul"=> (-Inf, Inf, -Inf, Inf),
                    "Wood"=> (4e3, 5e7, 1e-5, 4e-2),
                    "Churchill"=> (-Inf, Inf, -Inf, Inf)
)


function factor_friccion(Re,eD;función::Function=Colebrook)
    value = métodos["$función"]
    if value[1] <= Re <= value[2] && value[3] <= eD <= value[4]
        return función(Re,eD)
    else
        @warn "El método para el coeficiente de fricción seleccionado ($función), 
        solo es valido para $(value[1]) <= Re <= $(value[2]) y para $(value[3]) <= e/D <= $(value[4])"
        return función(Re,eD)
    end
end

function factor_friccion(Re;función::Function=laminar)
    value = métodos["$función"]
    if value[1] <= Re <= value[2]
        return función(Re)
    else
        @warn "El método para el coeficiente de fricción seleccionado ($función), 
        solo es valido para $(value[1]) <= Re <= $(value[2])"
        return función(Re)
    end
end    


end



