module MultivariateCreativeTelescoping

using StaticArrays
using Combinatorics
using Nemo
using DataStructures
using ExportAll

using FLINT_jll: libflint

global ord_ctr = 0; # used to create new types with metaprogramming
include("globalstats.jl")
include("primes.jl")
include("DataStructure_missing.jl")
# include("SLP.jl")

include("TypeCoefficient.jl")
include("TypeOreMonomial.jl")
include("TypeOrePolynomial.jl")
include("TypeMonOrder.jl")
include("TypeOreAlgebra.jl")
include("IndicialEqs.jl")

include("interface_FLINT.jl")

include("ParseInputOutput.jl")
include("OrePolyAddMul.jl")

include("Geobucket.jl")
include("Buchberger.jl")

include("symbolicpp.jl")
include("elimination.jl")
include("F4.jl")


include("TypeSigPair.jl")
include("F5.jl")

include("WeylClosure.jl")

include("echelon_form.jl")
include("Reductions.jl")
include("RepresentativeInIntMod.jl")
include("DerivRedMap.jl")
include("MCT.jl")

include("AlgebraMorphisms.jl")
include("CRT.jl")
include("CauchyInterpolation.jl")

include("annfs.jl")
include("testfunctions.jl")
# include("IndicialEqs.jl")

# include("multivariate_rational_interpolation.jl")
# include("multivariate_rational_interpolation2.jl")
include("multivariate_rational_interpolation3.jl")

@exportAll()

end
