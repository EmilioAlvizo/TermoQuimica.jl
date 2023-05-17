var documenterSearchIndex = {"docs":
[{"location":"api/#Api","page":"API","title":"Api","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"CurrentModule = TermoQuimica","category":"page"},{"location":"api/","page":"API","title":"API","text":"general.P_Antoine\ngeneral.T_Antoine","category":"page"},{"location":"api/#TermoQuimica.general.P_Antoine","page":"API","title":"TermoQuimica.general.P_Antoine","text":"P_Antoine(CA, T, i)\n\nlog_10P(mmHg)=A-fracBC+T(C)\n\nLos campos de entrada son:\n\nCA :: n×3 Array{Float,2} Constantes de Antoine para los n componentes\nT :: Float Es la temperatura del sistema en °C\ni :: Int Componente i al que se quiere conocer su presión\n\nSalida:\n\nRegresa la presión P en mmHg\n\nP :: m-element Array{Float64,1} presión\n\n\n\n\n\n","category":"function"},{"location":"api/#TermoQuimica.general.T_Antoine","page":"API","title":"TermoQuimica.general.T_Antoine","text":"T_Antoine(CA, P, i)\n\nT(C)=-C-fracBlog_10P(mmHg)-A\n\nLos campos de entrada son:\n\nCA :: n×3 Array{Float,2} Constantes de Antoine para los n componentes\nP :: Float Es la presión del sistema en mmHg\ni :: Int Componente i al que se quiere conocer su temperatura\n\nSalida:\n\nRegresa la temperatura T en °C\n\nT :: m-element Array{Float64,1} temperatura\n\n\n\n\n\n","category":"function"},{"location":"api/","page":"API","title":"API","text":"ideal.pxy\nideal.txy\nmargules.pxy\nmargules.txy\nvanlaar.pxy\nvanlaar.txy\nnrtl.pxy\nnrtl.txy\nwilson.pxy\nwilson.txy\nwilsonmod.pxy\nwilsonmod.txy\nuniquac.pxy\nuniquac.txy\nunifac.pxy\nunifac.txy\nunifacmod.pxy\nunifacmod.txy","category":"page"},{"location":"api/#TermoQuimica.ideal.pxy","page":"API","title":"TermoQuimica.ideal.pxy","text":"pxy(n,T,CA,xx;uni=u\"Torr\")\n\nLos campos de entrada son:\n\nn :: Int Es el numero de componentes\nT :: Float Es la temperatura del sistema\nCA :: n×3 Array{Float,2} Constantes de Antoine para los n componentes\nxx :: StepRange ó Vector Son los puntos en el liquido donde se buscara el equilibrio con el vapor\nuni Son las unidades en las que se desea el resultado, por defecto esta en Torr\n\nSalida:\n\nRegresa las fracciones de liquido x, las fracciones del vapor y y la presión P en Torr (m es la longitud de xx).\n\nx :: n×m Array{Float64,2} liquido\ny :: n×m Array{Float64,2} vapor\nP :: m-element Array{Float64,1} presión\n\n\n\n\n\n","category":"function"},{"location":"api/#TermoQuimica.ideal.txy","page":"API","title":"TermoQuimica.ideal.txy","text":"txy(n,Pobj,CA,xx;uni=u\"K\")\n\nLos campos de entrada son:\n\nn :: Int Es el numero de componentes\nPobj :: Float Es la presión del sistema\nCA :: n×3 Array{Float,2} Constantes de Antoine para los n componentes\nxx :: StepRange ó Vector Son los puntos en el liquido donde se buscara el equilibrio con el vapor\nuni Son las unidades en las que se desea el resultado, por defecto esta en grados kelvin\n\nSalida:\n\nRegresa las fracciones de liquido x, las fracciones del vapor y y la temperatura T en grados kelvin (m es la longitud de xx).\n\nx :: n×m Array{Float64,2} liquido\ny :: n×m Array{Float64,2} vapor\nT :: m-element Array{Float64,1} temperatura\n\n\n\n\n\n","category":"function"},{"location":"api/#TermoQuimica.margules.pxy","page":"API","title":"TermoQuimica.margules.pxy","text":"pxy(n,T,CA,Λ,xx;uni=u\"Torr\")\n\nLos campos de entrada son:\n\nn :: Int Es el numero de componentes\nT :: Float Es la temperatura del sistema\nCA :: n×3 Array{Float,2} Constantes de Antoine para los n componentes\nΛ :: n×n Array{Float,2} Parámetros de interacción entre los componentes en grados kelvin\nxx :: StepRange ó Vector Son los puntos en el liquido donde se buscara el equilibrio con el vapor\nuni Son las unidades en las que se desea el resultado, por defecto esta en Torr\n\nSalida:\n\nRegresa las fracciones de liquido x, las fracciones del vapor y y la presión P en Torr (m es la longitud de xx).\n\nx :: n×m Array{Float64,2} liquido\ny :: n×m Array{Float64,2} vapor\nP :: m-element Array{Float64,1} presión\n\n\n\n\n\n","category":"function"},{"location":"api/#TermoQuimica.margules.txy","page":"API","title":"TermoQuimica.margules.txy","text":"txy(n,Pobj,CA,Λ,xx;uni=u\"K\")\n\nLos campos de entrada son:\n\nn :: Int Es el numero de componentes\nPobj :: Float Es la presión del sistema\nCA :: n×3 Array{Float,2} Constantes de Antoine para los n componentes\nΛ :: n×n Array{Float,2} Parámetros de interacción entre los componentes en grados kelvin\nxx :: StepRange ó Vector Son los puntos en el liquido donde se buscara el equilibrio con el vapor\nuni Son las unidades en las que se desea el resultado, por defecto esta en grados kelvin\n\nSalida:\n\nRegresa las fracciones de liquido x, las fracciones del vapor y y la temperatura T en grados kelvin (m es la longitud de xx).\n\nx :: n×m Array{Float64,2} liquido\ny :: n×m Array{Float64,2} vapor\nT :: m-element Array{Float64,1} temperatura\n\n\n\n\n\n","category":"function"},{"location":"api/#TermoQuimica.vanlaar.pxy","page":"API","title":"TermoQuimica.vanlaar.pxy","text":"pxy(n,T,CA,Λ,xx;uni=u\"Torr\")\n\nLos campos de entrada son:\n\nn :: Int Es el numero de componentes\nT :: Float Es la temperatura del sistema\nCA :: n×3 Array{Float,2} Constantes de Antoine para los n componentes\nΛ :: n×n Array{Float,2} Parametros de interacción entre los componentes en grados kelvin\nxx :: StepRange ó Vector Son los puntos en el liquido donde se buscara el equilibrio con el vapor\nuni Son las unidades en las que se desea el resultado, por defecto esta en Torr\n\nSalida:\n\nRegresa las fracciones de liquido x, las fracciones del vapor y y la presión P en Torr (m es la longitud de xx).\n\nx :: n×m Array{Float64,2} liquido\ny :: n×m Array{Float64,2} vapor\nP :: m-element Array{Float64,1} presión\n\n\n\n\n\n","category":"function"},{"location":"api/#TermoQuimica.vanlaar.txy","page":"API","title":"TermoQuimica.vanlaar.txy","text":"txy(n,Pobj,CA,Λ,xx;uni=u\"K\")\n\nLos campos de entrada son:\n\nn :: Int Es el numero de componentes\nPobj :: Float Es la presión del sistema\nCA :: n×3 Array{Float,2} Constantes de Antoine para los n componentes\nΛ :: n×n Array{Float,2} Parametros de interacción entre los componentes en grados kelvin\nxx :: StepRange ó Vector Son los puntos en el liquido donde se buscara el equilibrio con el vapor\nuni Son las unidades en las que se desea el resultado, por defecto esta en grados kelvin\n\nSalida:\n\nRegresa las fracciones de liquido x, las fracciones del vapor y y la temperatura T en grados kelvin (m es la longitud de xx).\n\nx :: n×m Array{Float64,2} liquido\ny :: n×m Array{Float64,2} vapor\nT :: m-element Array{Float64,1} temperatura\n\n\n\n\n\n","category":"function"},{"location":"api/#TermoQuimica.nrtl.pxy","page":"API","title":"TermoQuimica.nrtl.pxy","text":"pxy(n,T,CA,Δg,α,xx;uni=u\"Torr\")\n\nLos campos de entrada son:\n\nn :: Int Es el numero de componentes\nT :: Float Es la temperatura del sistema\nCA :: n×3 Array{Float,2} Constantes de Antoine para los n componentes\nΔg :: n×n Array{Float,2} Parametros de interacción entre los componentes en grados kelvin\nα :: n×n Array{Float,2} Parametro de no aleatoriedad\nxx :: StepRange ó Vector Son los puntos en el liquido donde se buscara el equilibrio con el vapor\nuni Son las unidades en las que se desea el resultado, por defecto esta en Torr\n\nSalida:\n\nRegresa las fracciones de liquido x, las fracciones del vapor y y la presión P en Torr (m es la longitud de xx).\n\nx :: n×m Array{Float64,2} liquido\ny :: n×m Array{Float64,2} vapor\nP :: m-element Array{Float64,1} presión\n\n\n\n\n\n","category":"function"},{"location":"api/#TermoQuimica.nrtl.txy","page":"API","title":"TermoQuimica.nrtl.txy","text":"txy(n,Pobj,CA,Δg,α,xx;uni=u\"K\")\n\nLos campos de entrada son:\n\nn :: Int Es el numero de componentes\nPobj :: Float Es la presión del sistema\nCA :: n×3 Array{Float,2} Constantes de Antoine para los n componentes\nΔg :: n×n Array{Float,2} Parametros de interacción entre los componentes en grados kelvin\nα :: n×n Array{Float,2} Parametro de no aleatoriedad\nxx :: StepRange ó Vector Son los puntos en el liquido donde se buscara el equilibrio con el vapor\nuni Son las unidades en las que se desea el resultado, por defecto esta en grados kelvin\n\nSalida:\n\nRegresa las fracciones de liquido x, las fracciones del vapor y y la temperatura T en grados kelvin (m es la longitud de xx).\n\nx :: n×m Array{Float64,2} liquido\ny :: n×m Array{Float64,2} vapor\nT :: m-element Array{Float64,1} temperatura\n\n\n\n\n\n","category":"function"},{"location":"api/#TermoQuimica.wilson.pxy","page":"API","title":"TermoQuimica.wilson.pxy","text":"pxy(n,T,CA,λ,v,xx;uni=u\"Torr\")\n\nLos campos de entrada son:\n\nn :: Int Es el numero de componentes\nT :: Float Es la temperatura del sistema\nCA :: n×3 Array{Float,2} Constantes de Antoine para los n componentes\nλ :: n×n Array{Float,2} Parametros de interacción entre los componentes en grados kelvin\nv :: n-element Array{Float,1} Volumen molar\nxx :: StepRange ó Vector Son los puntos en el liquido donde se buscara el equilibrio con el vapor\nuni Son las unidades en las que se desea el resultado, por defecto esta en Torr\n\nSalida:\n\nRegresa las fracciones de liquido x, las fracciones del vapor y y la presión P en Torr (m es la longitud de xx).\n\nx :: n×m Array{Float64,2} liquido\ny :: n×m Array{Float64,2} vapor\nP :: m-element Array{Float64,1} presión\n\n\n\n\n\n","category":"function"},{"location":"api/#TermoQuimica.wilson.txy","page":"API","title":"TermoQuimica.wilson.txy","text":"txy(n,Pobj,CA,λ,v,xx;uni=u\"K\")\n\nLos campos de entrada son:\n\nn :: Int Es el numero de componentes\nPobj :: Float Es la presión del sistema\nCA :: n×3 Array{Float,2} Constantes de Antoine para los n componentes\nλ :: n×n Array{Float,2} Parametros de interacción entre los componentes en grados kelvin\nv :: n-element Array{Float,1} Volumen molar\nxx :: StepRange ó Vector Son los puntos en el liquido donde se buscara el equilibrio con el vapor\nuni Son las unidades en las que se desea el resultado, por defecto esta en grados kelvin\n\nSalida:\n\nRegresa las fracciones de liquido x, las fracciones del vapor y y la temperatura T en grados kelvin (m es la longitud de xx).\n\nx :: n×m Array{Float64,2} liquido\ny :: n×m Array{Float64,2} vapor\nT :: m-element Array{Float64,1} temperatura\n\n\n\n\n\n","category":"function"},{"location":"api/#TermoQuimica.wilsonmod.pxy","page":"API","title":"TermoQuimica.wilsonmod.pxy","text":"pxy(n,T,CA,a,b,c,v,xx;uni=u\"Torr\")\n\nLos campos de entrada son:\n\nn :: Int Es el numero de componentes\nT :: Float Es la temperatura del sistema\nCA :: n×3 Array{Float,2} Constantes de Antoine para los n componentes\na :: n×n Array{Float,2} Parametros de interacción a en fraccalmol\nb :: n×n Array{Float,2} Parametros de interacción b en fraccalmol K\nc :: n×n Array{Float,2} Parametros de interacción c en fraccalmol K^2\nv :: n-element Array{Float,1} Volumen molar\nxx :: StepRange ó Vector Son los puntos en el liquido donde se buscara el equilibrio con el vapor\nuni Son las unidades en las que se desea el resultado, por defecto esta en Torr\n\nSalida:\n\nRegresa las fracciones de liquido x, las fracciones del vapor y y la presión P en Torr (m es la longitud de xx).\n\nx :: n×m Array{Float64,2} liquido\ny :: n×m Array{Float64,2} vapor\nP :: m-element Array{Float64,1} presión\n\n\n\n\n\n","category":"function"},{"location":"api/#TermoQuimica.wilsonmod.txy","page":"API","title":"TermoQuimica.wilsonmod.txy","text":"txy(n,Pobj,CA,a,b,c,v,xx;uni=u\"K\")\n\nLos campos de entrada son:\n\nn :: Int Es el numero de componentes\nPobj :: Float Es la presión del sistema\nCA :: n×3 Array{Float,2} Constantes de Antoine para los n componentes\na :: n×n Array{Float,2} Parametros de interacción a en fraccalmol\nb :: n×n Array{Float,2} Parametros de interacción b en fraccalmol K\nc :: n×n Array{Float,2} Parametros de interacción c en fraccalmol K^2\nv :: n-element Array{Float,1} Volumen molar\nxx :: StepRange ó Vector Son los puntos en el liquido donde se buscara el equilibrio con el vapor\nuni Son las unidades en las que se desea el resultado, por defecto esta en grados kelvin\n\nSalida:\n\nRegresa las fracciones de liquido x, las fracciones del vapor y y la temperatura T en grados kelvin (m es la longitud de xx).\n\nx :: n×m Array{Float64,2} liquido\ny :: n×m Array{Float64,2} vapor\nT :: m-element Array{Float64,1} temperatura\n\n\n\n\n\n","category":"function"},{"location":"api/#TermoQuimica.uniquac.pxy","page":"API","title":"TermoQuimica.uniquac.pxy","text":"pxy(n,T,CA,Δu,q,r,xx;uni=u\"Torr\")\n\nLos campos de entrada son:\n\nn :: Int Es el numero de componentes\nT :: Float Es la temperatura del sistema\nCA :: n×3 Array{Float,2} Constantes de Antoine para los n componentes\nΔu :: n×n Array{Float,2} Parametros de interacción entre los componentes en grados kelvin\nq :: n-element Array{Float64,1} área de Van der Waals\nr :: n-element Array{Float64,1} volumen de Van der Waals\nxx :: StepRange ó Vector Son los puntos en el liquido donde se buscara el equilibrio con el vapor\nuni Son las unidades en las que se desea el resultado, por defecto esta en Torr\n\nSalida:\n\nRegresa las fracciones de liquido x, las fracciones del vapor y y la presión P en Torr (m es la longitud de xx).\n\nx :: n×m Array{Float64,2} liquido\ny :: n×m Array{Float64,2} vapor\nP :: m-element Array{Float64,1} presión\n\n\n\n\n\n","category":"function"},{"location":"api/#TermoQuimica.uniquac.txy","page":"API","title":"TermoQuimica.uniquac.txy","text":"txy(n,Pobj,CA,Δu,q,r,xx;uni=u\"K\")\n\nLos campos de entrada son:\n\nn :: Int Es el numero de componentes\nPobj :: Float Es la presión del sistema\nCA :: n×3 Array{Float,2} Constantes de Antoine para los n componentes\nΔu :: n×n Array{Float,2} Parametros de interacción entre los componentes en grados kelvin\nq :: n-element Array{Float64,1} área de Van der Waals\nr :: n-element Array{Float64,1} volumen de Van der Waals\nxx :: StepRange ó Vector Son los puntos en el liquido donde se buscara el equilibrio con el vapor\nuni Son las unidades en las que se desea el resultado, por defecto esta en grados kelvin\n\nSalida:\n\nRegresa las fracciones de liquido x, las fracciones del vapor y y la temperatura T en grados kelvin (m es la longitud de xx).\n\nx :: n×m Array{Float64,2} liquido\ny :: n×m Array{Float64,2} vapor\nT :: m-element Array{Float64,1} temperatura\n\n\n\n\n\n","category":"function"},{"location":"api/#TermoQuimica.unifac.pxy","page":"API","title":"TermoQuimica.unifac.pxy","text":"pxy(n,T,CA,com,xx;uni=u\"Torr\")\n\nLos campos de entrada son:\n\nn :: Int Es el numero de componentes\nT :: Float Es la temperatura del sistema\nCA :: n×3 Array{Float,2} Constantes de Antoine para los n componentes\ncom ::\nxx :: StepRange ó Vector Son los puntos en el liquido donde se buscara el equilibrio con el vapor\nuni Son las unidades en las que se desea el resultado, por defecto esta en Torr\n\nSalida:\n\nRegresa las fracciones de liquido x, las fracciones del vapor y y la presión P en Torr (m es la longitud de xx).\n\nx :: n×m Array{Float64,2} liquido\ny :: n×m Array{Float64,2} vapor\nP :: m-element Array{Float64,1} presión\n\n\n\n\n\n","category":"function"},{"location":"api/#TermoQuimica.unifac.txy","page":"API","title":"TermoQuimica.unifac.txy","text":"txy(n,Pobj,CA,com,xx;uni=u\"K\")\n\nLos campos de entrada son:\n\nn :: Int Es el numero de componentes\nPobj :: Float Es la presión del sistema\nCA :: n×3 Array{Float,2} Constantes de Antoine para los n componentes\ncom ::\nxx :: StepRange ó Vector Son los puntos en el liquido donde se buscara el equilibrio con el vapor\nuni Son las unidades en las que se desea el resultado, por defecto esta en grados kelvin\n\nSalida:\n\nRegresa las fracciones de liquido x, las fracciones del vapor y y la temperatura T en grados kelvin (m es la longitud de xx).\n\nx :: n×m Array{Float64,2} liquido\ny :: n×m Array{Float64,2} vapor\nT :: m-element Array{Float64,1} temperatura\n\n\n\n\n\n","category":"function"},{"location":"api/#TermoQuimica.unifacmod.pxy","page":"API","title":"TermoQuimica.unifacmod.pxy","text":"pxy(n,T,CA,com,xx;uni=u\"Torr\")\n\nLos campos de entrada son:\n\nn :: Int Es el numero de componentes\nT :: Float Es la temperatura del sistema\nCA :: n×3 Array{Float,2} Constantes de Antoine para los n componentes\ncom ::\nxx :: StepRange ó Vector Son los puntos en el liquido donde se buscara el equilibrio con el vapor\nuni Son las unidades en las que se desea el resultado, por defecto esta en Torr\n\nSalida:\n\nRegresa las fracciones de liquido x, las fracciones del vapor y y la presión P en Torr (m es la longitud de xx).\n\nx :: n×m Array{Float64,2} liquido\ny :: n×m Array{Float64,2} vapor\nP :: m-element Array{Float64,1} presión\n\n\n\n\n\n","category":"function"},{"location":"api/#TermoQuimica.unifacmod.txy","page":"API","title":"TermoQuimica.unifacmod.txy","text":"txy(n,Pobj,CA,com,xx;uni=u\"K\")\n\nLos campos de entrada son:\n\nn :: Int Es el numero de componentes\nPobj :: Float Es la presión del sistema\nCA :: n×3 Array{Float,2} Constantes de Antoine para los n componentes\ncom ::\nxx :: StepRange ó Vector Son los puntos en el liquido donde se buscara el equilibrio con el vapor\nuni Son las unidades en las que se desea el resultado, por defecto esta en grados kelvin\n\nSalida:\n\nRegresa las fracciones de liquido x, las fracciones del vapor y y la temperatura T en grados kelvin (m es la longitud de xx).\n\nx :: n×m Array{Float64,2} liquido\ny :: n×m Array{Float64,2} vapor\nT :: m-element Array{Float64,1} temperatura\n\n\n\n\n\n","category":"function"},{"location":"api/","page":"API","title":"API","text":"EOS.Van_Der_Waals\nEOS.Redlich_Kwong\nEOS.Soave_Redlich_Kwong\nEOS.Peng_Robinson\nsaft.pxy\nsaft.Mix","category":"page"},{"location":"api/#TermoQuimica.EOS.Van_Der_Waals","page":"API","title":"TermoQuimica.EOS.Van_Der_Waals","text":"Van_Der_Waals(T,P,Tc,Pc)\n\nEcuación de la forma:\n\nP=fracRTν-b-fracaν^2\n\nDonde:\n\nbeginaligned\na=frac27(RT_c)^264P_c  b=fracRT_c8P_c\nendaligned\n\nfracU-U^igRT=fracaRTν\n\nfracG-G^igRT=z-1-fracaRTν-log_10 (z-B)\n\nLos campos de entrada son:\n\nT :: Float Es la temperatura del sistema\nP :: Float Es la presión del sistema\nTc :: Float Es la temperatura critica del sistema\nPc :: Float Es la presión critica del sistema\n𝜔 :: Float Es el factor acéntrico \n\nSalida:\n\nRegresa un conjunto de propiedades termodinámicas.\n\nℍ :: Float Entalpía residual \n𝕌 :: Float Energía interna residual\n𝕊 :: Float Entropía residual\n𝔾 :: Float Energía libre de Gibbs residual\n𝔸 :: Float Energía libre de helmholtz residual\nφ :: Float Coeficiente de fugacidad \nf :: Float Fugacidad\nz :: Float Factor de compresibilidad \nν :: Float Volumen molar\n\n\n\n\n\n","category":"function"},{"location":"api/#TermoQuimica.EOS.Redlich_Kwong","page":"API","title":"TermoQuimica.EOS.Redlich_Kwong","text":"Redlich_Kwong(T,P,Tc,Pc)\n\nEcuación de la forma:\n\nP=fracRTν-b-fraca(T^15 ν(ν+b))\n\nDonde:\n\nbeginaligned\na=frac042748 R^2 T_c^25P_c  b=frac008664 RT_cP_c\nendaligned\n\nfracU-U^igRT=frac3alog_10 (1+zB)2bR T^15\n\nfracG-G^igRT=z-1+fracalog_10 (1+zB)bR T^15-log_10 (z-B)\n\nLos campos de entrada son:\n\nT :: Float Es la temperatura del sistema\nP :: Float Es la presión del sistema\nTc :: Float Es la temperatura critica del sistema\nPc :: Float Es la presión critica del sistema\n𝜔 :: Float Es el factor acéntrico \n\nSalida:\n\nRegresa un conjunto de propiedades termodinámicas.\n\nℍ :: Float Entalpía residual \n𝕌 :: Float Energía interna residual\n𝕊 :: Float Entropía residual\n𝔾 :: Float Energía libre de Gibbs residual\n𝔸 :: Float Energía libre de helmholtz residual\nφ :: Float Coeficiente de fugacidad \nf :: Float Fugacidad\nz :: Float Factor de compresibilidad \nv :: Float Volumen molar\n\n\n\n\n\n","category":"function"},{"location":"api/#TermoQuimica.EOS.Soave_Redlich_Kwong","page":"API","title":"TermoQuimica.EOS.Soave_Redlich_Kwong","text":"Soave_Redlich_Kwong(T,P,Tc,Pc,𝜔)\n\nEcuación de la forma:\n\nP=fracRTV-b-fracaαV(V+b)\n\nDonde:\n\nbeginaligned\na=frac042748 R^2 T_c^2P_c  b=frac008664 RT_cP_c  α=(1 + κ(1 -T_r^05))^2\nendaligned\n\nκ=048+1574w-0176w^2\n\nfracU-U^igRT=-frac(aα+asqrtT_rακ)log_10 (1+Bz)bRT\n\nfracS-S^igR=log_10(z-B)-fracasqrtT_rακlog_10(1+Bz)bRT\n\nLos campos de entrada son:\n\nT :: Float Es la temperatura del sistema\nP :: Float Es la presión del sistema\nTc :: Float Es la temperatura critica del sistema\nPc :: Float Es la presión critica del sistema\n𝜔 :: Float Es el factor acéntrico \n\nSalida:\n\nRegresa un conjunto de propiedades termodinámicas.\n\nℍ :: Float Entalpía residual \n𝕌 :: Float Energía interna residual\n𝕊 :: Float Entropía residual\n𝔾 :: Float Energía libre de Gibbs residual\n𝔸 :: Float Energía libre de helmholtz residual\nφ :: Float Coeficiente de fugacidad \nf :: Float Fugacidad\nz :: Float Factor de compresibilidad \nv :: Float Volumen molar\n\n\n\n\n\n","category":"function"},{"location":"api/#TermoQuimica.EOS.Peng_Robinson","page":"API","title":"TermoQuimica.EOS.Peng_Robinson","text":"Peng_Robinson(T,P,Tc,Pc,𝜔)\n\nEcuación de la forma:\n\nP=fracRTV-b-fracaαV^2 +2bV+b^2\n\nDonde:\n\nbeginaligned\na=frac045724 R^2 T_c^2P_c  b=frac00778 RT_cP_c  α=(1 + κ(1 -T_r^05))^2\nendaligned\n\nκ=037464 + 154226*w -026992*w^2\n\nfracU-U^igRT=-(aα+asqrtT_rακ)*(coth^-1sqrt2+tanh^-1frac(-1+Bz)sqrt2bRTsqrt2)\n\nfracG-G^igR=z-1-log(z-B)-aα*(coth^-1sqrt2+tanh^-1frac(B-z)zsqrt2bRTsqrt2)\n\nLos campos de entrada son:\n\nT :: Float Es la temperatura del sistema\nP :: Float Es la presión del sistema\nTc :: Float Es la temperatura critica del sistema\nPc :: Float Es la presión critica del sistema\n𝜔 :: Float Es el factor acéntrico \n\nSalida:\n\nRegresa un conjunto de propiedades termodinámicas.\n\nℍ :: Float Entalpía residual \n𝕌 :: Float Energía interna residual\n𝕊 :: Float Entropía residual\n𝔾 :: Float Energía libre de Gibbs residual\n𝔸 :: Float Energía libre de helmholtz residual\nφ :: Float Coeficiente de fugacidad \nf :: Float Fugacidad\nz :: Float Factor de compresibilidad \nv :: Float Volumen molar\n\n\n\n\n\n","category":"function"},{"location":"api/#TermoQuimica.saft.pxy","page":"API","title":"TermoQuimica.saft.pxy","text":"pxy(T,P,mix::mix,xx;uni=u\"Torr\")\n\nLos campos de entrada son:\n\nT :: Float Es la temperatura del sistema\nP :: m-element Array{Float64,1} presión inicial\nmix Datos de entrada de la mezcla\nxx :: StepRange ó Vector Son los puntos en el liquido donde se buscara el equilibrio con el vapor\nuni Son las unidades en las que se desea el resultado, por defecto esta en Torr\n\nSalida:\n\nRegresa las fracciones de liquido x, las fracciones del vapor y y la presión P en Torr (m es la longitud de xx).\n\nx :: n×m Array{Float64,2} liquido\ny :: n×m Array{Float64,2} vapor\nP :: m-element Array{Float64,1} presión\n\n\n\n\n\n","category":"function"},{"location":"api/#TermoQuimica.saft.Mix","page":"API","title":"TermoQuimica.saft.Mix","text":"mutable struct Mix\n\ngfhgfh\n\n\n\n\n\n","category":"type"},{"location":"api/","page":"API","title":"API","text":"accesorios.entrada_recta","category":"page"},{"location":"servidor/#Servidor","page":"Servidor","title":"Servidor","text":"","category":"section"},{"location":"servidor/","page":"Servidor","title":"Servidor","text":"","category":"page"},{"location":"servidor/#Puesta-en-marcha-del-servidor","page":"Servidor","title":"Puesta en marcha del servidor","text":"","category":"section"},{"location":"servidor/","page":"Servidor","title":"Servidor","text":"Si es usted el que hostea la pagina y debe levantar el servidor Genie, debe realizar los siguientes pasos:","category":"page"},{"location":"servidor/","page":"Servidor","title":"Servidor","text":"Diríjase a la carpeta bin la cual estará dentro de la carpeta de la pagina web, sabrá que a llegado a ella porque encontrara los archivos repl, server ... etc.\nEjecute el siguiente comando, que es para iniciar la aplicación y salir de la consola sin cerrar la aplicación:\n screen -S genie\nEjecute el archvo server.bat kaka:\n sh server.bat","category":"page"},{"location":"userguide/#Guia-de-usuario","page":"Guia de usuario","title":"Guia de usuario","text":"","category":"section"},{"location":"userguide/","page":"Guia de usuario","title":"Guia de usuario","text":"","category":"page"},{"location":"userguide/","page":"Guia de usuario","title":"Guia de usuario","text":"TermoQuímica es una pagina de apoyo para ayudar a realizar distintos tipos de cálculos útiles en la ingeniería química; En esta sección de la documentación, usted aprenderá como usar TermoQuímica, y descubrirá todo lo que puede hacer con ella durante un corto, pero divertido tutorial.","category":"page"},{"location":"userguide/","page":"Guia de usuario","title":"Guia de usuario","text":"Ahora, trataremos de resolver el siguiente problema.","category":"page"},{"location":"userguide/","page":"Guia de usuario","title":"Guia de usuario","text":"Elabore el diagrama P-X-Y para la mezcla etanol(1)-agua(2) a 273.15 K, usando el modelo de actividad de Van Laar; Para este ejercicio, las coeficientes de Van Laar son: A12 = 1.6798 y A21 = 0.9227, las constantes de Antoine para el etanol son: A = 8.12875, B = 1660.8713, C = 238.131 y para el agua A = 8.05573, B = 1723.6425, C = 233.08","category":"page"},{"location":"userguide/","page":"Guia de usuario","title":"Guia de usuario","text":"Lo primero que debemos hacer es hacer click en la pestaña de métodos y seleccionar Liquido/Vapor, escoger nuestro modelo que en este caso es Van Laar, después de hacer esto, aparecerá la siguiente pantalla, en la cual nos pide el numero de componentes, nuestra mezcla tiene dos componentes asi que ponemos 2 y damos a enviar.","category":"page"},{"location":"userguide/","page":"Guia de usuario","title":"Guia de usuario","text":"(Image: nocomp)","category":"page"},{"location":"userguide/","page":"Guia de usuario","title":"Guia de usuario","text":"Aparecerá la siguiente pantalla, la cual requiere que introduzcamos los distintos parámetros del modelo, que para el caso de Van Laar son los coeficientes, las constantes de Antoine, el rango en el que queremos evaluar la composición del componente 1, la temperatura y el tipo de diagrama, en este caso P-X-Y.","category":"page"},{"location":"userguide/","page":"Guia de usuario","title":"Guia de usuario","text":"(Image: nocomp)","category":"page"},{"location":"userguide/","page":"Guia de usuario","title":"Guia de usuario","text":"Habra que rellenar todos los campos y enviarlos.","category":"page"},{"location":"userguide/","page":"Guia de usuario","title":"Guia de usuario","text":"nota: Nota\nTenga en cuenta que para seleccionar un rango apropiado la sintaxis es inicio:paso:fin donde inicio y fin debe estar entre cero y uno. Para un sistema multicomponente los rangos deben ser iguales para todas las componentes, por ejemplo para una mezcla ternaria estos podrían ser los rangoscomponente 1 -> 0:0.05:1componente 2 -> 0:0.05:1","category":"page"},{"location":"userguide/","page":"Guia de usuario","title":"Guia de usuario","text":"(Image: nocomp)","category":"page"},{"location":"userguide/","page":"Guia de usuario","title":"Guia de usuario","text":"Una vez enviados los datos, se desplegaran los resultados en forma gráfica y en forma tabular (la forma gráfica solo se muestra para sistemas de dos o tres componentes), como se muestra en las siguientes imágenes.","category":"page"},{"location":"userguide/","page":"Guia de usuario","title":"Guia de usuario","text":"(Image: nocomp) (Image: nocomp)","category":"page"},{"location":"userguide/","page":"Guia de usuario","title":"Guia de usuario","text":"Al ver los resultados podríamos querer un numero mayor de puntos a calcular, para hacerlo, solo damos atrás en el navegador y cambiamos el rango del componente 1 de 0:0.1:1 a 0:0.01:1.","category":"page"},{"location":"userguide/","page":"Guia de usuario","title":"Guia de usuario","text":"(Image: nocomp)","category":"page"},{"location":"userguide/","page":"Guia de usuario","title":"Guia de usuario","text":"Y listo, ya tendríamos una gráfica mas bonita 🎉.","category":"page"},{"location":"userguide/","page":"Guia de usuario","title":"Guia de usuario","text":"(Image: nocomp)","category":"page"},{"location":"userguide/","page":"Guia de usuario","title":"Guia de usuario","text":"No solo eso podemos hacer, si en algún momento necesitáramos información acerca de puntos específicos, lo único que tenemos que hacer es escribirlos separados por una coma, de acuerdo con la imagen siguiente. ","category":"page"},{"location":"userguide/","page":"Guia de usuario","title":"Guia de usuario","text":"(Image: nocomp) (Image: nocomp)","category":"page"},{"location":"background/#Antecedentes","page":"Background","title":"Antecedentes","text":"","category":"section"},{"location":"background/","page":"Background","title":"Background","text":"","category":"page"},{"location":"background/#La-naturaleza-del-equilibrio","page":"Background","title":"La naturaleza del equilibrio","text":"","category":"section"},{"location":"background/","page":"Background","title":"Background","text":"El equilibrio es una condición en la que no se producen cambios en las propiedades macroscópicas de un sistema aislado con el tiempo. En el equilibrio, todos los potenciales que pueden causar cambios están exactamente equilibrados, por lo que no existe una fuerza motriz para ningún cambio en el sistema. Un sistema aislado que consiste en fases líquidas y vaporosas en contacto íntimo alcanza finalmente un estado final en el que no existe ninguna tendencia a que se produzcan cambios en el sistema. La temperatura, la presión y la composición de las fases alcanzan los valores finales que después permanecen fijos. El sistema está en equilibrio. Sin embargo, a nivel microscópico, las condiciones no son estáticas. Las moléculas que componen una fase en un instante dado no son las mismas moléculas que más tarde ocupan la misma fase. Las moléculas pasan constantemente de una fase a la otra. Sin embargo, la velocidad media de paso de las moléculas es la misma en ambas direcciones, y no se produce una transferencia neta de material entre fases. En la práctica de la ingeniería, la suposición de equilibrio se justifica cuando conduce a resultados de precisión satisfactoria. Por ejemplo, en la caldera de una columna de destilación, se suele suponer un equilibrio entre las fases de vapor y de líquido. En el caso de las tasas de vaporización finitas, se trata de una aproximación, pero no introduce un error significativo en los cálculos de ingeniería.","category":"page"},{"location":"background/#Equilibrio-Liquido/Vapor","page":"Background","title":"Equilibrio Liquido/Vapor","text":"","category":"section"},{"location":"background/","page":"Background","title":"Background","text":"Debido tanto a su tamaño como a su importancia, las columnas de destilación son características clave de las plantas químicas de todo el mundo. Las columnas de destilación son normalmente las columnas muy altas, de forma cilíndrica, que se ven en las plantas químicas que realizan la destilación, cuyo funcionamiento unitario separa los componentes de una mezcla en función de su punto de ebullición. Como era de esperar, la información sobre cómo las mezclas se separan en fases de vapor y líquido se convierte en la pieza clave de la información asociada a la destilación. ","category":"page"},{"location":"background/#Ejemplo","page":"Background","title":"Ejemplo","text":"","category":"section"},{"location":"background/","page":"Background","title":"Background","text":"Cuando se hierve una olla de agua a 1 atm, hay H2O líquido en la olla y vapor de H2O que sale de la olla, todo esto ocurre a 100 °C. En lugar de agua hirviendo sola (una sola sustancia, H2O), supongamos que tenemos una mezcla líquida equimolar de agua y otro compuesto, digamos etanol, en la olla a 1 atm. ¿A qué temperatura herviría la mezcla? ¿Estaría a 100 °C (el punto de ebullición normal del agua)? ¿Sería a 78 °C (el punto de ebullición normal del etanol)? ¿Estaría alrededor de 89 °C, que es un promedio ponderado de los puntos de ebullición normales basados en la composición (nuestro enfoque \"intuitivo\")? Además, ¿la composición del vapor que sale de la olla sería también equimolar?","category":"page"},{"location":"background/","page":"Background","title":"Background","text":"El experimento nos dice que una mezcla líquida equimolar de agua y etanol hervirá alrededor de 80 °C a 1 atm, mientras que la composición de la fase de vapor no será equimolar (como la del líquido), sino que será del 65% para el etanol. Estos resultados van en contra de nuestro enfoque intuitivo, pero también proporcionan una enorme posibilidad. Con solo hervir una mezcla líquida en una composición, podemos obtener una composición diferente en la fase de vapor. Así que en el ejemplo anterior, nuestra composición de vapor será enriquecida en etanol, en efecto \"separando\" parte del etanol del agua. Esta es la base de una de las operaciones unitarias más importantes y visibles en la ingeniería química: la columna de destilación. Obviamente, hay una gran necesidad de entender por qué nuestro enfoque intuitivo no funciona, así como de desarrollar una forma de predecir y/o modelar cuáles serán estos resultados para cualquier mezcla de vapor y líquido de interés en la industria de procesos químicos.","category":"page"},{"location":"background/","page":"Background","title":"Background","text":"Para comprender los resultados experimentales descritos en el párrafo anterior, tenemos que volver a los fenómenos de evaporación y ebullición para considerar la naturaleza del proceso a nivel molecular. Consideremos un líquido puro y comprimido en un contenedor de volumen fijo. El recipiente está dentro de un baño de calor que controla la temperatura y la mantiene constante. El líquido llena todo el espacio del contenedor. Si se deja solo, ¿se evaporará el líquido?. No, porque es un sistema cerrado y no hay espacio para que el vapor se llene. Sin embargo, digamos que ese recipiente está conectado a un segundo recipiente por un divisor y el contenido del segundo recipiente es evacuado (Figura 10-1), por lo que la presión del segundo recipiente es de 0 atm. Cuando el divisor es removido, ¿qué pasa? El líquido comenzará a expandirse en el espacio evacuado y una porción del líquido comenzará a vaporizarse porque cualquier burbuja de vapor que intente formarse desde la superficie del líquido no se detendrá (no hay presión, por lo tanto no hay fuerza para detenerlas) sobre el líquido en ese punto. Tanto la expansión como la vaporización del líquido continuarán hasta el punto en que parte del vapor comience a condensarse. Se establecerá un equilibrio dinámico donde la tasa de vaporización será igual a la tasa de condensación. La presión ejercida por este vapor sobre la superficie de este líquido se equilibra ahora con la presión ejercida en la fase líquida. Esta presión se denomina presión de vapor, que es la única presión a la que el líquido y el vapor pueden existir en equilibrio, y será una función de la temperatura del sistema.","category":"page"},{"location":"background/","page":"Background","title":"Background","text":"(Image: sdsdsd)","category":"page"},{"location":"background/","page":"Background","title":"Background","text":"¿Y si hiciéramos un experimento diferente? La configuración es inicialmente la misma, pero ahora tienes una mezcla de dos líquidos en lugar de un compuesto puro en el primer recipiente. El primer recipiente está de nuevo conectado a un segundo recipiente (cuyo contenido es evacuado) por un divisor. Quitamos el divisor y ¿qué pasa? ¿Es lo mismo que antes? Sí y no. Sí, hay expansión y vaporización y, en última instancia, se establece un equilibrio dinámico como antes. Pero, ¿cuál es la presión final en el sistema? Dado que ahora hay dos líquidos, ¿la presión final es sólo la suma de las presiones de vapor de ambos líquidos? No, no lo es. La presión de vapor de una sustancia pura es la presión sobre la superficie del líquido puro. Sin embargo, cuando se tiene una mezcla, hay más de un líquido presente en la superficie. Por lo tanto, para una superficie fija (como en este problema), habrá menos de cada componente en la interfaz líquido-vapor que si estuvieran allí por sí mismos. Por lo tanto, la presión parcial que cada componente ejercerá será diferente que si fuera un componente puro. ¿Cuánto menos? Bueno, hemos visto antes que cuando queremos calcular una propiedad de la mezcla, necesitamos usar las fracciones molares para dar cuenta del porcentaje de cada material en una mezcla","category":"page"},{"location":"background/","page":"Background","title":"Background","text":"Así que la pregunta es, \"¿Es la presión sobre la mezcla de líquido simplemente una suma de las presiones de vapor de los componentes puros (ponderadas por sus fracciones molares), como esto?\"","category":"page"},{"location":"background/","page":"Background","title":"Background","text":"Esto refleja nuestro enfoque \"intuitivo\" del capítulo anterior. La respuesta que dimos allí es la misma que ahora. Si el sistema se comporta como una solución ideal, esto sería el caso. Si no, entonces (como descubrimos de otras propiedades termofísicas en el capítulo 9) tendremos que dar cuenta de esos casos no ideales también.","category":"page"},{"location":"background/","page":"Background","title":"Background","text":"[2].","category":"page"},{"location":"background/#refs","page":"Background","title":"References","text":"","category":"section"},{"location":"background/","page":"Background","title":"Background","text":"[1] W. Tucker, Validated Numerics: A Short Introduction to Rigorous Computations, Princeton University Press (2011).","category":"page"},{"location":"background/","page":"Background","title":"Background","text":"[2] A. Haro, Automatic differentiation methods in computational dynamical systems: Invariant manifolds and normal forms of vector fields at fixed points, preprint.","category":"page"},{"location":"#TermoQuimica.jl","page":"Inicio","title":"TermoQuimica.jl","text":"","category":"section"},{"location":"","page":"Inicio","title":"Inicio","text":"Es un paquete construido en Julia, pensado para calcular el equilibrio liquido-vapor de diversos componentes con distintos métodos. Si busca mas información acerca del código pude visitar mi sitio en GitHub https://github.com/EmilioAlvizo/TermoQuimica.jl.git","category":"page"},{"location":"","page":"Inicio","title":"Inicio","text":"","category":"page"},{"location":"#Autor","page":"Inicio","title":"Autor","text":"","category":"section"},{"location":"","page":"Inicio","title":"Inicio","text":"Emilio Alvizo Velázquez, División de Ciencias e Ingenierías, Universidad de Guanajuato.","category":"page"}]
}
