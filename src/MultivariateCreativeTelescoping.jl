module MultivariateCreativeTelescoping

using StaticArrays
using Combinatorics
using Nemo
using DataStructures
using ExportAll

using FLINT_jll: libflint

include("globalstats.jl")
include("primes.jl")
include("DataStructure_missing.jl")

include("TypeCoefficient.jl")
include("TypeOreMonomial.jl")
include("TypeOrePolynomial.jl")
include("TypeMonOrder.jl")
include("TypeOreAlgebra.jl")

include("interface_FLINT.jl")

include("ParseInputOutput.jl")
include("OrePolyAddMul.jl")

include("Geobucket.jl")
include("Buchberger.jl")

include("symbolicpp.jl")
include("elimination.jl")
include("fglm.jl")
include("F4.jl")


include("TypeSigPair.jl")
include("F5.jl")

include("WeylClosure.jl")

include("echelon_form.jl")
include("Reductions.jl")
include("RepresentativeInIntMod.jl")
include("DerivRedMap.jl")
include("DerivRedMapMany.jl")
include("MCT.jl")
include("MCTMany.jl")

include("AlgebraMorphisms.jl")
include("CRT.jl")
include("CauchyInterpolation.jl")

include("annfs.jl")
include("dfinite_parser/database_LDE.jl")
include("dfinite_parser/ann_sum_prod.jl")
include("dfinite_parser/minimal_polynomial.jl")
include("dfinite_parser/dfinite_algebraic_leaf.jl")
include("dfinite_parser/ann_comp_right_rat.jl")
include("dfinite_parser/ann_comp_right_alg.jl")
include("dfinite_parser/dfinite_hyperexp_division.jl")
include("dfinite_parser/ann_poly_power.jl")
include("dfinite_parser/dfinite_parser.jl")
include("testfunctions.jl")
# include("IndicialEqs.jl")

# include("multivariate_rational_interpolation.jl")
# include("multivariate_rational_interpolation2.jl")
include("multivariate_rational_interpolation3.jl")
include("multivariate_rational_interpolation_known_support.jl")

@exportAll()

end
