module TermoQuimica

using Unitful
using NLsolve

export EOS, general, vanlaar, ideal, margules, wilson, wilsonmod, nrtl, uniquac, unifac, unifacmod, saft
include("EOS.jl")
include("general.jl")
include("ideal.jl")
include("margules.jl")
include("vanlaar.jl")
include("wilson.jl")
include("wilsonmod.jl")
include("nrtl.jl")
include("uniquac.jl")
include("unifac.jl")
include("unifacmod.jl")
include("saft.jl")

end # module
