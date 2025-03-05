struct F4Param{A,B,C,D} <: GBParam end 

f4_param(;geobucket :: Val{A} = Val(false), 
         stophol :: Val{B} = Val(false),
         stat :: Val{C} = Val(false),
         debug :: Val{D} = Val(false)) where {A,B,C,D} = F4Param{A,B,C,D}() 

geobucket(:: F4Param{A,B,C,D}) where {A,B,C,D} = A 
stophol(:: F4Param{A,B,C,D}) where {A,B,C,D} = B 
stat(:: F4Param{A,B,C,D}) where {A,B,C,D} = C 
debug(::F4Param{A,B,C,D}) where {A,B,C,D} = D



# Building upon AbsOreMonomial define monomials indexed by integers

struct IdxMon  <: AbsOreMonomial{Int}
    ind :: Int
end


function mon(P :: OrePoly{T,M}, i :: Integer) where {T, M<:IdxMon}
    return mons(P)[i].ind 
end

Base.@propagate_inbounds Base.getindex(P :: OrePoly{T,M}, i :: Integer) where {T,M<:IdxMon}= (getindex(P.coeffs, i), getindex(P.mons, i).ind)
Base.@propagate_inbounds function Base.setindex!(P::OrePoly{T, M}, t::Tuple{T, I}, i) where {M<:IdxMon, T,I<:Int}
    setindex!(P.coeffs, t[1], i)
    setindex!(P.mons, IdxMon(t[2]), i)
end




Base.@propagate_inbounds function monomialsubs(P :: OrePoly, dic :: Dict{M,IdxMon}) where {M <: AbsOreMonomial}
    co = copy(coeffs(P))
    mo = Vector{IdxMon}(undef,length(P))
     for i in 1:length(P)
        mo[i] = dic[mon(P,i)]
    end
    return OrePoly(co,mo)
end

function monomialsubsinv(P :: OrePoly, v :: Vector{M}) where {M <: AbsOreMonomial}
    co = copy(coeffs(P))
    mo = [v[i.ind] for i in mons(P)]
    return OrePoly(co,mo)
end



#. PIVOT MATRIX

# The structure computed by symbolic preprocessing
mutable struct F4Matrix{M,I,T,alg <: OreAlg}
    A :: alg

    # Translation from column indices to monomial indices and conversely
    idxtomonomial::Vector{M}
    monomialtoidx::Dict{M,I}

    nbcolumns::Int
    # εcolumn::I           # first column corresponding to an ε term. beyond this
    #                      # column, we only carry the information but mostly
    #                      # ignore it.

    rows::Vector{OrePoly{T, I}} # sparse rows

    pivots::Vector{Int}         # pivots[j] is 0 or the index of a row whose
                                # first nonzero entry is in column j

    donotpivot::BitSet    # The rows that may not be used as pivots

    newpivots::BitSet

    inputrows::BitSet
end

Base.@propagate_inbounds function ispivot(mx::F4Matrix, i)
    row = mx.rows[i]
    return @inbounds !isempty(row) && mx.pivots[mon(row,1)] == i
end

function row(mx::F4Matrix, i :: Int)
    return monomialsubsinv(mx.rows[i], mx.idxtomonomial)
end


#Append `p` to `rows` and add all monomials of `p` in `todo`, expect the ε-terms and
#the monomials that are already in `done`


function pushrow!(A::alg,
                  rows::Vector{OrePoly{K,M}},
                  p::OrePoly{K,M},
                  todo::Set{M},
                  done::Set{M},
                  ) where {K,M, alg <: OreAlg}

    push!(rows, p)

    for k in 1:length(p)
        m = mon(p, k)
        if !(m ∈ done)
            # if isepsilon(ctx.mo, m)
            #     # we don't look for reducers for ε terms
            #     push!(done, m)
            # else
                push!(todo, m)
            # end
        end
    end
end





#For each monomial in `todo`, find a reducer in `basis` and add it to `rows` with `pushrow!`

function saturate!(A::alg,
                   rows::Vector{OrePoly{K,M}},
                   basis::Vector{OrePoly{K,M}},
                   todo::Set{M}, done::Set{M},
                   param ::F4Param,
                   geob :: GeoBucket) where {K,M, alg <: OreAlg}

    # recursively find all possible reducers

    while !isempty(todo)
        m = pop!(todo)
        m ∈ done && continue
        push!(done, m)

        # isepsilon(ctx.mo, m) && continue

        # search for a reducer
        # TODO benchmark and improve
        for red in Iterators.reverse(basis)
            lmred = mon(red,1)
            if divide(lmred, m,A)
                if geobucket(param)
                    sred = shift(red, m,A,geob)
                else 
                    sred = shift(red, m,A)
                end
                if stat(param) 
                    globalstats.counters[:f4_nb_reducer_computed] += 1
                    globalstats.counters[:f4_size_reducer] += length(sred)
                    globalstats.counters[:f4_size_m]+= sum(m/lmred)
                    globalstats.counters[:f4_deg_reducer] += maxdeg(sred)
                    add_reducer_globalstats!(i, m/lmred)
                end

                pushrow!(A, rows, sred, todo, done)
                break
            end
        end
    end

    return
end





#Construct the matrix from a list of sparse row.
#Concretely, it computes the column indices, reindex all rows so that the
#index is now the column index, not the index in the monomial table; it finds
#the pivots and normalize them if required.

function f4matrix(A::alg,
                  rows_::Vector{OrePoly{K,M}},
                  monomialset::Set{M}, donotpivot::BitSet,
                  inputrows::BitSet
                  ) where {K,M, alg <: OreAlg}

    monomials = collect(M, monomialset)
    sort!(monomials, order=order(A), rev=true)

    monomialtoidx = Dict{M,IdxMon}()
    for (i, m) in enumerate(monomials)
        monomialtoidx[m] = IdxMon(i)
    end

    # Better to do it in place?
    rows = Vector{OrePoly{K,IdxMon}}(undef,length(rows_))
    for i in 1:length(rows_)
        rows[i] = monomialsubs(rows_[i], monomialtoidx)
    end 
        # Among all possible pivots for a column, we choose the last
    pivots = zeros(Int, length(monomials))
    for (i, row) in enumerate(rows)
        #IdxMon(i) ∈ donotpivot && continue
        pivots[mon(row,1)] = i
    end

    for p in pivots
        if p != 0
            row = rows[p]
            lc = coeff(row, 1)
            if !isone(lc,ctx(A))
                # normalize if required
                mult = inv(lc, ctx(A))
                newco = [mul(mult, c,ctx(A)) for c in row.coeffs]
                rows[p] = OrePoly(newco,mons(row))
            end
        end
    end

    # εcolumn = findlast(m -> !isepsilon(ctx.mo, m), monomials) + 1

    f4matrix = F4Matrix{M, IdxMon,eltype_co(A), typeof(A)}(
        A, monomials, monomialtoidx,
        length(monomials), rows, pivots,
        donotpivot, BitSet(), inputrows)

end



#Input:
#    - `ctx`, polynomial context
#    - `ix`, context for indexing monomials
#    - `pols`, list of indexed polynomials
#    - `basis`, list of unindexed polynomials, to be used as reducers
#    - `interreduction`, if false, then input rows are never used as pivots
#    - `spairreduction`, if true, then no reducers are searched for the leading terms
#Output:
#    - A F4Matrix datastructure.

#NB: the reducers are appended to rows, and the input polynomials are reindexed.


function symbolicpp(A::alg,
    pols::Vector{OrePoly{K,M}},
    basis::Vector{OrePoly{K,M}},
    param ::F4Param,
    geob :: GeoBucket,
    interreduction::Bool=true,
    spairrecution::Bool=false,
   ) where {alg <: OreAlg, K, M}

    isempty(pols) && error("nothing to preprocess")

    rows = OrePoly{K,M}[]

    todo = Set{M}()
    done = Set{M}()
    inputrows = BitSet()
    for (i, p) in enumerate(pols)
        if spairrecution
            # in F4, we reduce together a bunch of half-S-pairs. In this case,
            # we do not need to search a reducer for leading terms, because it
            # is already reduced by the other half.
            push!(done, mon(p,1))
        end
        pushrow!(A, rows, p, todo, done)
        push!(inputrows, length(rows))
    end

    saturate!(A, rows, basis, todo, done,param,geob)

    return f4matrix(A, rows, done, interreduction ? Set{M}() : inputrows, inputrows)
end

function interreductionmx(A :: OreAlg,
                          basis::Vector{OrePoly{I, T}},
                          param :: F4Param,
                          geob :: GeoBucket
                          ) where {I, T}


    # increasing order of leading monomial

    basis = sort(basis, order=order(A))
    for p in basis 
        if length(p) == 0 
            error("zero vector in basis")
        end
    end
    rows = OrePoly{I, T}[]
    currentbasis = eltype(basis)[]

    todo = Set{T}()
    done = Set{T}()
    inputrows = BitSet()

    for (i, p) in enumerate(basis)
        pushrow!(A, rows, p, todo, done)
        push!(inputrows, length(rows))
        saturate!(A, rows, currentbasis, todo, done,param,geob)
        push!(currentbasis, p)
    end

    return f4matrix(A, rows, done, BitSet(), inputrows)
end


function shift(P :: OrePoly{K,M}, m :: M, A :: alg,geob :: GeoBucket) where {K,M,alg <: OreAlg}
    addmul_geobucket!(geob,one(ctx(A)),m/mon(P,1),P,A)
    return normalform(geob,A)
end
