module MultivariateCreativeTelescoping

using StaticArrays
using Combinatorics
using Nemo
using DataStructures
using ExportAll

include("globalstats.jl")
include("primes.jl")

include("TypeCoefficient.jl")
include("TypeOreMonomial.jl")
include("TypeOrePolynomial.jl")
include("TypeMonOrder.jl")
include("TypeOreAlgebra.jl")

include("ParseInputOutput.jl")
include("OrePolyAddMul.jl")

include("Buchberger.jl")

include("symbolicpp.jl")
include("elimination.jl")
include("F4.jl")

include("TypeSigPair.jl")
include("F5.jl")

include("WeylClosure.jl")

include("RepresentativeInIntMod.jl")
include("DerivRedMap.jl")
include("MCT.jl")

include("AlgebraMorphisms.jl")
include("CRT.jl")
include("CauchyInterpolation.jl")

include("testfunctions.jl")

@exportAll()

end
