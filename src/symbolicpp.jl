struct F4Param{A,B,C,D,E,F} <: GBParam end 

f4_param(;geobucket :: Val{A} = Val(false), 
         stophol :: Val{B} = Val(false),
         stat :: Val{C} = Val(false),
         debug :: Val{D} = Val(false),
         select_reducer :: Val{E} = Val(:last),
         tracer :: Val{F} = Val(:none)) where {A,B,C,D,E,F} = F4Param{A,B,C,D,E,F}() 

geobucket(:: F4Param{A,B,C,D,E,F}) where {A,B,C,D,E,F} = A 
stophol(:: F4Param{A,B,C,D,E,F}) where {A,B,C,D,E,F} = B 
stat(:: F4Param{A,B,C,D,E,F}) where {A,B,C,D,E,F} = C 
debug(::F4Param{A,B,C,D,E,F}) where {A,B,C,D,E,F} = D
select_reducer(::F4Param{A,B,C,D,E,F}) where {A,B,C,D,E,F} = E
# the two functions below check whether a tracer is used
learn(::F4Param{A,B,C,D,E,F}) where {A,B,C,D,E,F} = F == :learn
apply(::F4Param{A,B,C,D,E,F}) where {A,B,C,D,E,F} = F == :apply



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
            push!(todo, m)
        end
    end
end

function pushrow_delayed!(A::alg,
                  p::OrePoly{K,M},
                  todo::Set{M},
                  done::Set{M},
                  ) where {K,M, alg <: OreAlg}

    for k in 1:length(p)
        m = mon(p, k)
        if !(m ∈ done)
            push!(todo, m)
        end
    end
end





# For each monomial in `todo`, find a reducer in `basis` and add it to `rows` with `pushrow!`

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

        # search for a reducer

        i = select_reducer(A,basis,m,Val(select_reducer(param)))
        if i == 0 # no reducer found
            continue 
        end
        red = basis[i]
        if geobucket(param)
            sred = shift(red, m,A,geob)
        else 
            sred = shift(red, m,A)
        end
        if stat(param) 
            lmred = mon(red,1)
            globalstats.counters[:f4_nb_reducer_computed] += 1
            globalstats.counters[:f4_size_reducer] += length(sred)
            globalstats.counters[:f4_size_m]+= sum(m/lmred)
            globalstats.counters[:f4_deg_reducer] += maxdeg(sred)
            add_reducer_globalstats!(i, m/lmred)
        end

        pushrow!(A, rows, sred, todo, done)
    end

    return
end



function saturate2!(A::alg,
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
        for i in length(basis):-1:1
            red = basis[i]
            lmred = mon(red,1)
            if divide(lmred, m,A)
                @assert select_reducer(A,basis,m,Val(select_reducer(param))) == i 
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

function select_reducer(A :: OreAlg,
               basis :: Vector{OrePoly{K,M}},
               m :: OreMonVE,
               :: Val{E}) where {K,M,E}
    we= 2^31
    j = 0 
    for i in 1:length(basis)
        if divide(mon(basis[i],1),m,A) 
            w = weight(basis[i],Val(E))
            if w <= we
                j = i 
                we = w 
            end
        end
    end
    return j
end

function weight(g :: OrePoly, ::Val{:last}) 
    return 2^31
end

function weight(g :: OrePoly, :: Val{:shortest})
    return length(g)
end

function weight(g :: OrePoly, :: Val{:elim})
    w = 0 
    d = sum(mon(g,1)) 
    for m in mons(g) 
        s = sum(m)
        if d > s 
            w += d - s 
        else 
            w +=1
        end
    end
    return w 
end

function weight(g :: OrePoly, :: Val{:elim_coeffsize})
    w = 0
    d = sum(mon(g,1)) 
    for m in mons(g) 
        s = sum(m)
        if d > s 
            w += d -s 
        else 
            w += 1
        end
    end
    return w*weight(denominator(g)) 
end

function weight(g :: OrePoly, :: Val{:coeffsize})
    return weight(denominator(g)) 
end


weight(c :: UInt32) = 32 
weight(c :: UnivRatFunModp) = 32*length(c)
weight(c :: UnivRatFunQQ) = weight(numerator(c)) + weight(Nemo.denominator(c))
weight(c :: RatFunModp) = 32*length(c)
weight(c :: RatFunQQ) = weight(numerator(c)) + weight(Nemo.denominator(c))
weight(c :: ZZMPolyRingElem) = sum(weight(cc) for cc in coefficients(c))
weight(c :: ZZPolyRingElem) = sum(weight(cc) for cc in coeffs(c))
weight(c :: ZZRingElem) = ceil(Int,log2(abs(c)))
weight(c :: QQFieldElem) = weight(numerator(c)) + weight(Nemo.denominator(c))


function div2!(f :: OrePoly, g :: Vector{OrePoly{K,M}} ,A :: OreAlg; full :: Val{B} = Val(true),param :: GBParam = DefaultParam()) where {K,M,B}
    r=1 
    while r <= length(f) 
        m = mon(f,r)
        i = select_reducer(A,g,m,Val(:shortest))
        if i > 0 
            f = reduce!(f,r, g[i], A,param = param)
        else 
            r = r + 1
        end
    end
    return f
end




# Construct the matrix from a list of sparse row.
# Concretely, it computes the column indices, reindex all rows so that monomials
# are replaced by its associated column index; it finds
# the pivots and normalize them if required.

function f4matrix(A::alg,
                  rows_::Vector{OrePoly{K,M}},
                  monomialset::Set{M}, donotpivot::BitSet,
                  inputrows::BitSet, newpivots :: BitSet
                  ) where {K,M, alg <: OreAlg}

    # sorts monomials and create a bijection with a finite subset of the integers
    monomials = collect(M, monomialset)
    sort!(monomials, order=order(A), rev=true)

    monomialtoidx = Dict{M,IdxMon}()
    for (i, m) in enumerate(monomials)
        monomialtoidx[m] = IdxMon(i)
    end

    # replace monomials by their associated column index
    rows = Vector{OrePoly{K,IdxMon}}(undef,length(rows_))
    for i in 1:length(rows_)
        rows[i] = monomialsubs(rows_[i], monomialtoidx)
    end 

    # Among all possible pivots for a column, we choose the last, except for newpivots who take priority
    
    pivots = zeros(Int, length(monomials))
    for (i, row) in enumerate(rows)
        IdxMon(i) ∈ donotpivot && continue
        pivots[mon(row,1)] = i
    end
    for i in newpivots
        row = rows[i] 
        pivots[mon(row,1)] = i 
    end

    # normalize pivots
    for p in pivots
        if p != 0
            row = rows[p]
            lc = coeff(row, 1)
            if !isone(lc,ctx(A))
                makemonic!(row,A)
            end
        end
    end

    f4matrix = F4Matrix{M, IdxMon,eltype_co(A), typeof(A)}(
        A, monomials, monomialtoidx,
        length(monomials), rows, pivots,
        donotpivot, newpivots, inputrows)

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
    spairrecution::Bool=false
   ) where {alg <: OreAlg, K, M}

    isempty(pols) && error("nothing to preprocess")

    rows = OrePoly{K,M}[]

    todo = Set{M}()
    done = Set{M}()
    inputrows = BitSet()
    newpivots = BitSet()
    for (i, p) in enumerate(pols)
        if spairrecution
            if stat(param)
                globalstats.counters[:f4_nb_reducer_computed] += 1
            end
            # in F4, we reduce together a bunch of half-S-pairs. In this case,
            # we do not need to search a reducer for leading terms, because it
            # is already reduced by the other half.
            push!(done, mon(p,1))
        else 
            # otherwise the spair could be a new pivot
            m = mon(p,1)
            if !any(divide(mon(g,1),m,A) for g in basis) && !any(divide(mon(g,1),m,A) for g in rows)
                push!(done, m)
                push!(newpivots,length(rows)+1)
            end
            pushrow!(A, rows, p, todo, done)
        end
        push!(inputrows, length(rows))
    end

    saturate!(A, rows, basis, todo, done,param,geob)

    return f4matrix(A, rows, done, interreduction ? BitSet() : inputrows, inputrows, newpivots)
end

function interreductionmx(A :: OreAlg,
                          basis::Vector{OrePoly{I, T}},
                          param :: F4Param,
                          geob :: GeoBucket
                          ) where {I, T}


    # sort by increasing order of leading monomials
    basis = sort(basis, order=order(A))

    rows = OrePoly{I, T}[]
    currentbasis = eltype(basis)[]

    todo = Set{T}()
    done = Set{T}()
    inputrows = BitSet()

    for p in basis
        pushrow!(A, rows, p, todo, done)
        push!(inputrows, length(rows))
        saturate!(A, rows, currentbasis, todo, done,param,geob)
        push!(currentbasis, p)
    end

    return f4matrix(A, rows, done, BitSet(), inputrows, BitSet())
end


function shift(P :: OrePoly{K,M}, m :: M, A :: alg,geob :: GeoBucket) where {K,M,alg <: OreAlg}
    addmul_geobucket!(geob,one(ctx(A)),m/mon(P,1),P,A)
    return normalform(geob,A)
end
