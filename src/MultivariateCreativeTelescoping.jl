module MultivariateCreativeTelescoping

using StaticArrays
using Combinatorics
using Nemo
#using DataStructures
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


@exportAll()

end
