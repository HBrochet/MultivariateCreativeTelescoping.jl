
mutable struct Spair{M}
    left :: Int                 # index of the left part in the current basis
    right :: Int                # right part
    lcm :: M                    # lcm of the leading monomials of both parts
    lcmisprod :: Bool           # is the lcm the product of the leading monomials and does l and r commute ? 
    redundant :: Bool           # a flag used in Gebauer-Moller algorithm
end

# # S-pairs are ordered by lcm
function Base.Order.lt(ord :: O, a :: Spair{T}, b :: Spair{T}, ) where {T,O <: AbsMonomialOrder}
    lt(ord, a.lcm, b.lcm)
end

struct F4Trace 
    spairs :: Vector{Vector{Tuple{Int64, Int64}}} # spairs[i] is the list of spairs to consider during the ith loop of f4 
end
f4_trace() = F4Trace(Vector{Tuple{Int64,Int64}}[])



# All data that we need during the F4 algorithm
mutable struct PartialGB{T,M, Alg,B,C,D}
    # must be indexed monomials (ie monomials are indices)

    alg :: Alg 
    basis :: Vector{OrePoly{T,M}}          # current basis
    active :: BitSet            # indices of active elements in basis
    newrels :: Vector{OrePoly{T,M}}        # new relations for which the S-pairs are not yet generated
    candidates :: Vector{OrePoly{T,M}}     # halfpairs to be reduced
    spairs :: Vector{Spair{M}} # spairs
    param :: F4Param{B,C,D}
    geob :: GeoBucket{T, M}
    function PartialGB{T, M, Alg,B,C,D}(
        A :: Alg, generators :: Vector{OrePoly{T,M}},param :: F4Param{B,C,D}
    ) where {T, M, Alg <: OreAlg,B,C,D}
        gens = copy(generators)
        filter!(p -> !isempty(p), gens)
        geob = GeoBucket(zero(A))
        gens2 = interreduce(A,gens,param,geob)
        new(A, OrePoly{T,M}[],
            BitSet(), gens2, OrePoly{T,M}[], Vector{Spair{M}}(),param,geob)
    end
end


function spair(pgb :: PartialGB, i :: Int, j :: Int,)
    p = pgb.basis[i]
    q = pgb.basis[j]

    thelcm = lcm(mon(p,1), mon(q,1))

    # product rule ignored as it is almost never useful

    # if (thelcm == mon(p,1)*mon(q,1)) &&  commute(p,q,pgb.alg)
    #     prod = true 
    # else 
    #     prod = false
    # end
    prod = false

    Spair{eltype_mo(pgb.alg)}(i, j, thelcm, prod, false)
end


function partialgb(generators, A :: OreAlg,param :: F4Param{B,C,D}) where{B,C,D}
    PartialGB{eltype_co(A), eltype_mo(A), typeof(A),B,C,D}(A, generators,param)
end


iscomplete(pgb :: PartialGB) = isempty(pgb.newrels) && isempty(pgb.spairs) && isempty(pgb.candidates)

function pushrel!(pgb :: PartialGB, rel :: OrePoly) 
    # Implements Gebaur-Möller criterion for removing redundant Spairs
    # See Mora 2005, Solving Pol. Eq. Sys. II, Fig. 25.1
    # See also Lemma 25.1.9

    push!(pgb.basis, rel)
    debug(pgb.param) && @debug "New relation with leading monomial $(mon(rel,1).exp)"
    s = length(pgb.basis)
    if isempty(pgb.active)
        push!(pgb.active, s)
        return
    end

    M = eltype_mo(pgb.alg)
    SpairT = Spair{M}

    lmrel = mon(rel,1)
    newsp = Dict{Int,SpairT}()
    for i in 1:s-1
        if iscompatible(lmrel, mon(pgb.basis[i],1),pgb.alg)
            newsp[i] = spair(pgb, i, s)
        end
    end

    for sp in pgb.spairs
        sp.redundant && continue
        if divide(lmrel, sp.lcm,pgb.alg) && newsp[sp.left].lcm != sp.lcm && newsp[sp.right].lcm != sp.lcm
            sp.redundant = true
        end
    end

    for sp in values(newsp)
        if !(sp.left ∈ pgb.active)
            sp.redundant = true
        end
    end

    for spi in values(newsp)
        spi.redundant && continue
        for spj in values(newsp)
            spj.redundant && continue
            if spi.lcm != spj.lcm && divide(spj.lcm, spi.lcm, pgb.alg)
                spi.redundant = true
            end
        end
    end

    spbylcm = Dict{M, SpairT}()
    dontconsider = Set{M}()
    ctr = 0 # number of pair eliminated with product criterion
    for sp in values(newsp)
        sp.redundant && continue
        if sp.lcmisprod
            ctr += 1
            push!(dontconsider, sp.lcm)
            delete!(spbylcm, sp.lcm)
        elseif !(sp.lcm ∈ dontconsider)
            # Among all S-pairs with identical cm, select those with smallest left number.
            prevsp = get(spbylcm, sp.lcm, sp)
            if prevsp.left >= sp.left
                spbylcm[sp.lcm] = sp
            end
        end
    end

    for sp in values(spbylcm)
        push!(pgb.spairs, sp)
    end
    if stat(pgb.param)
        globalstats.counters[:f4_candidate_spairs] += length(spbylcm)
        globalstats.counters[:f4_eliminated_spairs_with_GM] += s-1 - length(spbylcm) - ctr 
        globalstats.counters[:f4_eliminated_spairs_with_prod_crit] += ctr
    end


    for i in pgb.active
        # otherwise, lmrel is not really reduced
        # @assert !divides(mγ, leadingmonomial(pgb.basis[i]), lmrel)
        # Actually this can happen because pgb.basis[i] may be also a new relation, so the current new relation
        # may have not been reduced by multiples of pgb.basis[i]
        if divide(lmrel, mon(pgb.basis[i],1),pgb.alg)
            delete!(pgb.active, i)
        end
    end

    push!(pgb.active, s)
    return
end


function generatespairs!(pgb :: PartialGB{M,T}) where {M,T}
    debug(pgb.param) && @debug "generating s-pairs for $(length(pgb.newrels)) new relations"

    for rel in pgb.newrels
        isempty(rel) && continue
        makemonic!(rel,pgb.alg)
        pushrel!(pgb, rel)
    end
    empty!(pgb.newrels)
    return
end

function add_new_rels!(pgb :: PartialGB{M,T}) where {M,T}
    debug(pgb.param) && @debug "adding $(length(pgb.newrels)) new relations"

    for rel in pgb.newrels
        isempty(rel) && continue
        debug(pgb.param) && @debug "New relation with leading monomial $(mon(rel,1).exp)"

        makemonic!(rel,pgb.alg)

        push!(pgb.basis, rel)
        lmrel = mon(rel,1)
        for i in pgb.active
            if divide(lmrel, mon(pgb.basis[i],1),pgb.alg)
                delete!(pgb.active, i)
            end
        end
        push!(pgb.active, length(pgb.basis))
    end
    empty!(pgb.newrels)
    return 
end

function selectspairs!(pgb :: PartialGB, param :: F4Param)
    isempty(pgb.spairs) && return

    if learn(param) 
        spairs = Tuple{Int,Int}[]
    end


    #todo: is it the best strategy for an elim order ? 
    sort!(pgb.spairs, rev=true, order=order(pgb.alg))
    deg = max_deg_block(last(pgb.spairs).lcm,pgb.alg)

    ctr = 0 
    while !isempty(pgb.spairs) && ctr < 1000
        sp = last(pgb.spairs)
        if max_deg_block(sp.lcm,pgb.alg) > deg 
            break 
        end
        l = pgb.basis[sp.left]
        r = pgb.basis[sp.right]
        lcm = sp.lcm 
        A = pgb.alg
        if geobucket(param)
            geob = pgb.geob
            addmul_geobucket!(geob,lc(r),lcm/lm(l),l,A) 
            addmul_geobucket!(geob,opp(lc(l),ctx(A)),lcm/lm(r),r,A) 
            p = normalform(geob, A)
        else
            p = sub(mul(lcm/lm(l),l,A),
                    mul(lcm/lm(r),r,A),
                    A)
        end
        if !isempty(p)
            push!(pgb.candidates, p)
        end

        if learn(param) 
            push!(spairs, (sp.left, sp.right))
        end
        pop!(pgb.spairs)
        ctr += 1 
    end
    if learn(param) 
        return spairs
    else
        return 
    end
end

function selectspairs!(pgb :: PartialGB, param :: F4Param, trace ::F4Trace, ctr :: Int)
    for (i,j) in trace.spairs[ctr]
        l = pgb.basis[i]
        r = pgb.basis[j]
        thelcm = lcm(mon(l,1), mon(r,1))
        A = pgb.alg
        if geobucket(param)
            geob = pgb.geob
            addmul_geobucket!(geob,lc(r),thelcm/lm(l),l,A) 
            addmul_geobucket!(geob,opp(lc(l),ctx(A)),thelcm/lm(r),r,A) 
            p = normalform(geob, A)
        else
            p = sub(mul(thelcm/lm(l),l,A),
                    mul(thelcm/lm(r),r,A),
                    A)
        end
        if !isempty(p)
            push!(pgb.candidates, p)
        end
    end
    return 
end



#
function selectspol!(pgb :: PartialGB)
    isempty(pgb.spairs) && return zero(pgb.alg)

    sort!(pgb.spairs, rev=true, order=order(pgb.alg))
    
    sp =  pop!(pgb.spairs)

    l = pgb.basis[sp.left]
    r = pgb.basis[sp.right]
    lcm = sp.lcm 
    if stat(pgb.param) 
        globalstats.counters[:bbg_nb_spair] += 1
    end
    res = sub(mul(lcm/lm(l),l,pgb.alg),
             mul(lcm/lm(r),r,pgb.alg),
             pgb.alg)
    return res
end


function findnewrels!(pgb :: PartialGB)
    isempty(pgb.candidates) && return
   
    activebasis = [pgb.basis[i] for i in pgb.active]

    rows = pgb.candidates

    debug(pgb.param) && @debug "symbolic preprocessing"
    mx = symbolicpp(pgb.alg, rows, activebasis, pgb.param,pgb.geob,false, false)

    debug(pgb.param) && @debug "reducing a sparse $(length(mx.rows))x$(mx.nbcolumns) matrix"
    reducepivots!(mx,pgb.param)
    reduce!(mx,pgb.param)
    # todo are the new pivots nicely reduced ? is it not already done
    reducenewpivots!(mx,pgb.param)
    pgb.newrels = [row(mx, r) for r in mx.newpivots]
    empty!(pgb.candidates)
    return
end

function findnewrels!(pgb :: PartialGB, spairs :: Vector{Tuple{Int,Int}})
    isempty(pgb.candidates) && return
    non_zero_sps = Tuple{Int,Int}[] # remember spairs that were not reduced to zero
    activebasis = [pgb.basis[i] for i in pgb.active]

    rows = pgb.candidates

    debug(pgb.param) && @debug "symbolic preprocessing"
    mx = symbolicpp(pgb.alg, rows, activebasis, pgb.param,pgb.geob,false, false)

    debug(pgb.param) && @debug "reducing a sparse $(length(mx.rows))x$(mx.nbcolumns) matrix"
    reducepivots!(mx,pgb.param)
    reduce!(mx,pgb.param)
    reducenewpivots!(mx,pgb.param)
    
    # fill non_zero_sps
    for i in mx.inputrows
        if ispivot(mx,i) 
            push!(non_zero_sps, spairs[i])
        end
    end
    pgb.newrels = [row(mx, r) for r in mx.newpivots]
    empty!(pgb.candidates)
    return non_zero_sps
end

"""
    f4(gens :: Vector{OrePoly{T,M}}, A :: Alg)

Return a reduced Gröbner basis of the left ideal generated by gens for the monomial order defined in A.
"""
function f4(gens :: Vector{OrePoly{K,M}}, A :: Alg; param :: F4Param = f4_param(), trace ::F4Trace = f4_trace()) where {K,M,Alg <: OreAlg}
    isempty(gens) && return gens
    pgb = partialgb(gens, A,param)
    round = 0
    if apply(param) && length(trace.spairs) == 0
        sort!(basis, lt = (x,y) -> lt(order(A),x[1][2], y[1][2]), rev = true)
        return pgb.basis
    end

    while !iscomplete(pgb)
        round += 1
        if stophol(param) && isholonomic(pgb.basis,A)
            debug(param) && @debug "The current partial gb generates a holonomic ideal"
            delete_spairs_with_T!(pgb,A)
        end

        debug(param) && @debug "starting round $round, there are $(length(pgb.spairs)) spairs to be treated"
        
        if learn(param)
            generatespairs!(pgb)
            spairs = selectspairs!(pgb,param)
            non_zero_spairs = findnewrels!(pgb,spairs)
            push!(trace.spairs, non_zero_spairs)
            @assert length(non_zero_spairs) == length(pgb.newrels)
        elseif apply(param)                 
            add_new_rels!(pgb)
            round > length(trace.spairs) && break
            nb_spairs = length(trace.spairs[round])
            selectspairs!(pgb,param,trace,round)
            findnewrels!(pgb)
            @assert nb_spairs == length(pgb.newrels) 
        else 
            generatespairs!(pgb)
            selectspairs!(pgb,param)
            findnewrels!(pgb)
        end

        debug(param) && @debug "found $(length(pgb.newrels)) new relations"    
    end

    basis = [pgb.basis[i] for i in pgb.active]
    stophol(param) && isholonomic(basis,A) && delete_op_with_T!(basis,A)
    basis = interreduce(A, basis,param,pgb.geob)
    # reducebasis!(basis,A)
    sort!(basis, lt = (x,y) -> lt(order(A),x[1][2], y[1][2]), rev = true)
    if learn(param) 
        filter!(p -> length(p) > 0,trace.spairs)
        return basis, trace
    else
        return basis
    end
end

# wrapper for using crt with tracer = Val(:apply) 
function f4(trace ::F4Trace, gens :: Vector{OrePoly{K,M}}, A :: Alg, param :: F4Param) where {K,M,Alg <: OreAlg}
    par = f4param_learn_to_apply(param)
    return f4(gens,A, param= par, trace = trace)
end
# wrapper for using crt with tracer = Val(:learn) 

function f4(gens :: Vector{OrePoly{K,M}}, A :: Alg,param :: F4Param; tracer :: Val{B} = Val(false)) where {B,K,M,Alg <: OreAlg}
    B && @assert learn(param) # for compatibility with compute_with_crt
    return f4(gens,A, param= param)
end



function f4_mri(A :: OreAlg,gen :: Vector{OrePoly{K,M}}, param :: F4Param) where {K,M}
    # flatten the result to allow rational reconstruction with mri 
    tmp =  f4(gen,A,param = param)
    res = tmp[1]
    for i in 2:length(tmp)
        append!(res,tmp[i])
    end
    return res 
end




function Buchberger2(gens_ :: Vector{OrePoly{K,M}}, A :: Alg;param ::F4Param = f4_param()) where {K,M,Alg <: OreAlg}
    gens = deepcopy(gens_)
    pgb = partialgb(gens, A,param)
    round = 0

    while !iscomplete(pgb)
        if stophol(param) && isholonomic(pgb.basis,A)
            debug(param) && @debug "The current partial gb generates a holonomic ideal"
            delete_spairs_with_T!(pgb,A)
        end
        round += 1
        debug(param) && @debug "starting round $round, there are $(length(pgb.spairs)) spairs to be treated"

        generatespairs!(pgb)
        spol = selectspol!(pgb)
        if iszero(spol) 
            debug(param) && @debug "S-polynom is zero"
            continue
        end
        
        debug(param) && @debug "reducing S-polynom with lm $(lm(spol).exp)"
        spol = div2!(spol,pgb.basis,A,param=param)

        if !iszero(spol)
            makemonic!(spol,A)
            push!(pgb.newrels,spol) 
        end


        debug(param) && @debug "found $(length(pgb.newrels)) new relation"
    end

    # globalstats.timings[:interreduction] += @elapsed begin
    #     basis = interreduce(ctx, [pgb.basis[i] for i in pgb.active])
    # end
    # basis = interreduce(A, [pgb.basis[i] for i in pgb.active])
    basis = [pgb.basis[i] for i in pgb.active] # todo: check what pgb.active is 
    stophol(param) && isholonomic(basis,A) && delete_op_with_T!(basis,A)
    reducebasis!(basis,A)
    sort!(basis, lt = (x,y) -> lt(order(A),x[1][2], y[1][2]), rev = true)
    return basis
end

function Buchberger2_mri(A :: OreAlg,gen :: Vector{OrePoly{K,M}}, param :: F4Param) where {K,M}
    # flatten the result to allow rational reconstruction with mri 
    tmp =  Buchberger2(gen,A,param = param)
    res = tmp[1]
    for i in 2:length(tmp)
        append!(res,tmp[i])
    end
    return res 
end



function delete_spairs_with_T!(pgb :: PartialGB, A :: OreAlg)
    N = nvars(A)
    hasT = Int[]
    for (i,g) in enumerate(pgb.basis) 
        if mon(g,1)[N] > 0 
            push!(hasT,i)
        end
    end
    todelete = Int[] 
    for (i,s) in enumerate(pgb.spairs)
        if (s.left in hasT) || (s.right in hasT)
            push!(todelete,i)
        end
    end

    reverse!(todelete)
    for s in todelete
        deleteat!(pgb.spairs,s)
    end
end

function delete_op_with_T!(v :: Vector{OrePoly{C, M}},A:: OreAlg) where {C,M}
    N = nvars(A)
    todelete = Int[]
    for (i,g) in enumerate(v) 
        if mon(g,1)[N] > 0 
            push!(todelete,i)
        end
    end
    reverse!(todelete)
    for s in todelete
        deleteat!(v,s)
    end
end




function print_pgb(pgb :: PartialGB)
    println("basis")
    prettyprint(pgb.basis,pgb.alg)
    println("actice")
    println(pgb.active)
    println("newrels")
    prettyprint(pgb.newrels,pgb.alg)
    println("candidates")
    prettyprint(pgb.candidates,pgb.alg)
    println("spairs")
    println(pgb.spairs)
end