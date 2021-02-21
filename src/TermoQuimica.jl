module TermoQuimica

using Unitful
using NLsolve

export EOS, general, vanlaar, ideal, margules, wilson, wilsonmod, nrtl, uniquac, unifac, unifacmod
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


#=
tolueno 7.13657  1457.2871  231.827
acetone 7.31742  1315.6735  240.479
water   8.05573  1723.6425  233.08 

uniquac
i j     aij       aji 
1 2     115.14      29.585
1 3     814.64      334.88
2 3     317.21      -22.287

r1 3.9228
r2 2.5735
r3 0.92
q1 2.968
q2 2.336
q3 1.4

=#

end # module
