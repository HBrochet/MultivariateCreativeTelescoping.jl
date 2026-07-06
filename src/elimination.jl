include("elimination_field.jl")
include("elimination_ring.jl")

# #. Interreduction of Gröbner bases

function interreduce(alg :: OreAlg,
                      pols :: Vector{OrePoly{I, T}},
                      param :: F4Param,
                      geob :: GeoBucket
                      ) where {I, T}

    mx = interreductionmx(alg, pols,param,geob)

    reducepivots!(mx,param)
    reduce!(mx,param)

    pols = [row(mx, i) for i in mx.inputrows]

    filter!(p -> !isempty(p), pols)
    sort!(pols, order=order(alg))

    return pols
end

