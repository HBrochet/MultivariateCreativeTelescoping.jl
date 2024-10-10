function init_gb_sat(g :: Vector{OrePoly{K,M}},bound :: Int, i :: Int, A :: OreAlg) where {K,M}
    res = copy(g)
    res2 = OrePoly{K,M}[]


    T = makepoly(one(ctx(A)),makemon(i,A))

    for (j,p) in enumerate(g)
        mulT = mul(T,p,A)

        while mon(mulT,1)[i] <= bound 
            push!(res,mulT)
            mulT = mul(T,mulT,A)
        end
        push!(res2, mulT)
    end
    return res,res2 
end

function update_gb_sat!(gb :: Vector{OrePoly{K,M}}, lastpow :: Vector{OrePoly{K,M}},bound :: Int, i :: Int, A :: OreAlg) where {K,M}
    T = makepoly(one(ctx(A)),makemon(i,A))
    for j in 1:length(lastpow)
        while mon(lastpow[j],1)[i] <= bound 
            push!(gb,lastpow[j])
            lastpow[j] = mul(T,lastpow[j],A)
        end
    end
end


# i is the index of the saturation variable in A
function saturation(g_ :: Vector{OrePoly{T,M}}, i :: Int, A :: OreAlg; bnd::Integer = 1,morestep::Integer=0,stophol::Bool=true,f4b::Bool=false) where {T,M}
    bound = bnd # bound on the degree in the loc variable
    g = deepcopy(g_)

    tmp = mul(A.pols_loc[1],parse_OrePoly(A.inp.locvars[1][1],A),A)
    tmp = sub(tmp,one(A),A)
    push!(g, tmp)
    gb, lastpow = init_gb_sat(f5(g,A), bound,i,A)
    if f4b 
        gb = f4(gb,A)
    else 
        gb = f5(gb,A)
    end
    res = filter(p -> mon(p,1)[i] == 0, gb) 
    @debug "gb initialized, it contains $(length(gb)) vectors"
    hol = isholonomic(res,A)
    while !hol || morestep > 0 
        bound = bound + 1
        if hol 
            morestep -= 1 
        end
        update_gb_sat!(gb,lastpow,bound,i,A)
        if f4b 
            gb = f4(gb,A;stophol=stophol)
        else 
            gb = f5(gb,A,stophol=stophol)
        end
        res = filter(p -> mon(p,1)[i] == 0, gb) 
        hol = isholonomic(res,A)
    end
    return res 
end


struct SingLocInit{T,C,M,O}
    A :: OreAlg{T,C,M,O}
end

function singular_locus_init(A :: OreAlg)
    return SingLocInit(to_ore_alg_with_rat_coeffs(A))
end
Base.show(io::IO, SLI :: SingLocInit) = print(io,"SingLocInit")

"""
Computes a polynomial vanishing on the singular locus of a D-module presented by a set of relation gens in an algebra A
"""
function singular_locus(gens :: Vector{OrePoly{T,M}}, A :: OreAlg,init :: SingLocInit) where {T,M}
    nA = init.A
    ngens = map_algebras(gens,A,nA)
    gb = f5(ngens,nA)
    clear_denominators!(gb,nA)

    pol = ctx(nA).R(1)
    for p in gb 
        pol = lcm(pol, numerator(coeff(p,1)))
    end
    return pol
end

struct WeylClosureInit{T,C,M,T2,C2,M2,O1,O2}
    SLI :: SingLocInit{T,C,M,O1}
    A :: OreAlg{T2,C2,M2,O2}
end
Base.show(io::IO, :: WeylClosureInit) = print(io,"WeylClosureInit")

"""
    weyl_closure_init(A :: OreAlg)

Initialisation function for weyl_closure
"""
function weyl_closure_init(A :: OreAlg)
    return WeylClosureInit(singular_locus_init(A),add_loc_var("_T","1",A))
end

function weyl_closure_internal(gens::Vector{OrePoly{T,M}}, A :: OreAlg, init :: WeylClosureInit, stophol :: Bool) where {T,M}
    nA = init.A
    sl = singular_locus(gens,A,init.SLI)
    if sl == one(A)
        return f5(gens,A)
    else 
        nemo_p = prod(fact[1] for fact in factor_squarefree(sl) )
        p = OrePoly([nemo_p],[makemon(-1,A)])
    end
    # change algebra so that _T corresponds to 1/p
    (polsloc,diff_pols_loc) = compute_locpol_and_derivs([mystring(p,A)],nA)
    empty!(nA.pols_loc)
    append!(nA.pols_loc,polsloc)
    empty!(nA.diff_pols_loc) 
    append!(nA.diff_pols_loc, diff_pols_loc) 
    ngens = map_algebras(gens,A,nA)
    gb = saturation(ngens,nA.nrdv+2*nA.npdv +nA.npv+nA.nlv,nA; stophol=stophol)
    return map_algebras(gb,nA,A)
end

"""
    weyl_closure(gens::Vector{OrePoly{T,M}}, A :: OreAlg, init :: WeylClosureInit)

Return a Gröbner basis of a submodule of the Weyl closure of the module generated by gens.
"""
function weyl_closure(gens::Vector{OrePoly{T,M}}, A :: OreAlg, init :: WeylClosureInit; stophol ::Bool =true) where {T,M}
    if ctx(A) isa UnivRatFunModpCtx || ctx(A) isa RatFunModpCtx || ctx(A) isa UnivRatFunModPCtx || ctx(A) isa RatFunModPCtx || ctx(A) isa NmodΓ
        return weyl_closure_internal(gens,A, init, stophol)
    else 
        return compute_with_CRT(weyl_closure_internal,A, gens, A, init, stophol)
    end
end