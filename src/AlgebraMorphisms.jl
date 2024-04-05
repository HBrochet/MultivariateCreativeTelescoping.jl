"""
returns a new Ore Algebra where every polynomial variables associated to differential variables
are converted to rational variables.
"""
function to_ore_alg_with_rat_coeffs(A :: OreAlg)
    inp = deepcopy(A.inp)
    # remove polynomial variables in order 
    for s in inp.poldiffvars[1]
        inp.order = replace(inp.order," "*s*" "=>" "; count=1)
    end
    # change variables 
    append!(inp.ratdiffvars[1],inp.poldiffvars[1])
    append!(inp.ratdiffvars[2],inp.poldiffvars[2])
    empty!(inp.poldiffvars[1])
    empty!(inp.poldiffvars[2])
    return OreAlg(inp)
end


# """
# Convert an operator p in the algebra A to a new operator in the algebra nA returned by the function 
# to_ore_alg_with_rat_coeffs.

# """
# function pol_coeffs_to_rat_coeffs(p :: OrePoly, A :: OreAlg, nA :: OreAlg)
#     s = mystring(p,A)
#     return parse_OrePoly(s,nA)
# end

""" 
Returns a new algebra built from A by adding one localisation variable s associated to the polynomial p.
the variable s is also added to nomul.
"""
function add_loc_var(s :: String, p :: String, A ::OreAlg)
    inp = deepcopy(A.inp)
    inp.order = "lex "*s*" > "*inp.order
    push!(inp.locvars[1],s)
    push!(inp.locvars[2],p)
    push!(inp.nomul,s)
    return OreAlg(inp)
end

# """
# Convert an operator p in the algebra A to a new operator in the algebra nA returned by the function 
# add_loc_var.

# """
# function morph_add_loc_var(p :: OrePoly, A:: OreAlg, nA :: OreAlg)

"""
Given a polynomial p in the algebra A, returns its representation in the new algebra nA. 
This function assume that p can be represented in nA
"""
function map_algebras(p :: OrePoly, A:: OreAlg,nA :: OreAlg)
    s = mystring(p,A)
    return parse_OrePoly(s,nA)
end

function map_algebras(gens :: Vector{OrePoly{T,M}},A:: OreAlg,nA :: OreAlg) where {T,M}
    return [map_algebras(p,A,nA) for p in gens]
end



""" change the characteristic of the algebra"""
function change_alg_char(prime :: Int, A :: OreAlg)
    if typeof(ctx(A)) isa RatFunQQ
        return change_alg_char_ratfun(prime, A)
    else
        error("not implemented yet")
    end
end


function change_alg_char_ratfun(prime_ :: Integer, A :: OreAlg)
    prime = Int(prime_)
    #Context 
    CRing = Native.GF(prime)
    R = parent(change_coefficient_ring(CRing,ctx(A).vars[1]))
    F = fraction_field(R)
    nvars = [change_coefficient_ring(CRing,var) for var in ctx(A).vars]
    nctx = RatFunModpCtx(F,R,nvars,prime)

    ev_ratvars = Dict{String,RatFunModp}()
    for k in keys(A.ratvars)
        ev_ratvars[k] = F(change_coefficient_ring(CRing, numerator(A.ratvars[k])))
    end
    inp = deepcopy(A.inp)
    inp.char = prime

    T = RatFunModp
    M = eltype_mo(A)
    tmpA = OreAlg{eltype1_ctx(nctx),typeof(nctx),eltype_mo(A),}(A.strvar_to_indexp,
                                                                A.indexp_to_strvar,
                                                                ev_ratvars,                                                               
                                                                A.rvars_to_int,
                                                                A.int_to_rvars,
                                                                A.nrdv,
                                                                A.npdv,
                                                                A.npv,
                                                                A.nlv,
                                                                OrePoly{T,M}[],
                                                                Vector{OrePoly{T,M}}[],
                                                                A.nomul,
                                                                A.order,
                                                                nctx,
                                                                A.inp)

    pl = [change_coefficient_field(CRing,p,tmpA) for p in A.pols_loc]
    dpl = [change_coefficient_field(CRing,v,tmpA) for v in A.diff_pols_loc]

    return OreAlg{eltype1_ctx(nctx),typeof(nctx),eltype_mo(A),}(A.strvar_to_indexp,
                                                                A.indexp_to_strvar,
                                                                ev_ratvars,                                                               
                                                                A.rvars_to_int,
                                                                A.int_to_rvars,
                                                                A.nrdv,
                                                                A.npdv,
                                                                A.npv,
                                                                A.nlv,
                                                                pl,
                                                                dpl,
                                                                A.nomul,
                                                                A.order,
                                                                nctx,
                                                                A.inp)
end

function change_alg_char_QQ(prime :: Integer, A :: OreAlg)
    #Context 
    nctx = Nmod32Γ(prime)

    nratvars = Dict{String,UInt32}()

    inp = deepcopy(A.inp)
    inp.char = prime

    T = UInt32
    M = eltype_mo(A)
    tmpA = OreAlg{eltype1_ctx(nctx),typeof(nctx),eltype_mo(A),}(A.strvar_to_indexp,
                                                                A.indexp_to_strvar,
                                                                nratvars,                                                               
                                                                A.rvars_to_int,
                                                                A.int_to_rvars,
                                                                A.nrdv,
                                                                A.npdv,
                                                                A.npv,
                                                                A.nlv,
                                                                OrePoly{T,M}[],
                                                                Vector{OrePoly{T,M}}[],
                                                                A.nomul,
                                                                A.order,
                                                                nctx,
                                                                A.inp)

    pl = [change_coefficient_field(p,tmpA) for p in A.pols_loc]
    dpl = [change_coefficient_field(v,tmpA) for v in A.diff_pols_loc]

    return OreAlg{eltype1_ctx(nctx),typeof(nctx),eltype_mo(A),}(A.strvar_to_indexp,
                                                                A.indexp_to_strvar,
                                                                nratvars,                                                               
                                                                A.rvars_to_int,
                                                                A.int_to_rvars,
                                                                A.nrdv,
                                                                A.npdv,
                                                                A.npv,
                                                                A.nlv,
                                                                pl,
                                                                dpl,
                                                                A.nomul,
                                                                A.order,
                                                                nctx,
                                                                A.inp)
end

"""change the coefficient field of an OrePoly from the fraction field of a Nemo ring to the fraction field of a new Nemo ring"""
function change_coefficient_field(R :: Ring, gs :: Vector{OrePoly{T,M}}, A :: OreAlg) where {T, M}
    return [change_coefficient_field(R,g,A) for g in gs]
end


function change_coefficient_field(R :: Ring, p :: OrePoly, A :: OreAlg)
    cs = [change_coefficient_field(R,c,A) for c in coeffs(p)]
    return OrePoly(cs,deepcopy(mons(p)))
end

function change_coefficient_field(R :: Ring, p :: Generic.FracFieldElem{ZZMPolyRingElem}, A :: OreAlg)
    num = numerator(p)
    den = Nemo.denominator(p)
    F = ctx(A).F
    return F(change_coefficient_ring(R,num)) / F(change_coefficient_ring(R,den))
end

function change_coefficient_field(:: Ring,a :: Any,:: OreAlg)
    return a 
end

function change_coefficient_field(R :: Ring,wci :: WeylClosureInit, A :: OreAlg)
    return WeylClosureInit(change_coefficient_field(R,wci.SLI,A), change_coefficient_field(R,wci.A,A))
end

function change_coefficient_field(R :: Ring,sli :: SingLocInit, A :: OreAlg)
    return SingLocInit(change_coefficient_field(R,sli.A,A))
end


function change_coefficient_field(R :: Ring,A :: OreAlg, args...)
    return (change_coefficient_field(R,arg,A) for arg in args)
end

function change_coefficient_field(:: Ring, vA ::OreAlg, A ::OreAlg)
    return change_alg_char_ratfun(char(A),vA)
end


# change coefficient field from QQ to modp 
function change_coefficient_field(A :: OreAlg, args...)
    return (change_coefficient_field(arg,A) for arg in args)
end

function change_coefficient_field(vA ::OreAlg, A ::OreAlg)
    return change_alg_char_QQ(char(A),vA)
end

function change_coefficient_field(sli :: SingLocInit, A :: OreAlg)
    return SingLocInit(change_alg_char_ratfun(char(A),sli.A)
    )
end

function change_coefficient_field(wci :: WeylClosureInit, A :: OreAlg)
    return WeylClosureInit(change_coefficient_field(wci.SLI,A), change_coefficient_field(wci.A,A))
end

function change_coefficient_field(gs :: Vector{OrePoly{T,M}}, A :: OreAlg) where {T, M}
    return [change_coefficient_field(g,A) for g in gs]
end


function change_coefficient_field(p :: OrePoly, A :: OreAlg)
    cs = [change_coefficient_field(c,A) for c in coeffs(p)]
    return OrePoly(cs,deepcopy(mons(p)))
end

function change_coefficient_field(p :: QQFieldElem, A :: OreAlg)
    num = p.num
    den = p.den
    return mul(convertn(num,ctx(A)), inv(convertn(den,ctx(A)), ctx(A)), ctx(A))
end

function change_coefficient_field(a :: Any,:: OreAlg)
    return a 
end

"""
Given as input an algebra A  with only one parameter and a point p, return a new algebra where the parameter was evaluated to p.
"""
function evaluate_parameter_algebra(p :: Int,A :: OreAlg)
    new_ctx = Nmod32Γ(ctx(A).char)
    T = eltype1_ctx(new_ctx)
    M = eltype_mo(A)
    tmpA = OreAlg{T,typeof(new_ctx),eltype_mo(A)}(A.strvar_to_indexp,
                    A.indexp_to_strvar,
                    Dict{String,UInt32}(),
                    Dict{String,Int}(),
                    String[],
                    A.nrdv,
                    A.npdv,
                    A.npv,
                    A.nlv,
                    OrePoly{T,M}[],
                    Vector{OrePoly{T,M}}[],
                    A.nomul,
                    order(A),
                    new_ctx,
                    A.inp)
    ev_pol_locs = evaluate_parameter(A.pols_loc, p, tmpA)
    ev_diff_pols_loc = evaluate_parameter(A.diff_pols_loc, p, tmpA)
    return OreAlg{eltype1_ctx(new_ctx),typeof(new_ctx),eltype_mo(A)}(A.strvar_to_indexp,
                    A.indexp_to_strvar,
                    Dict{String,UInt32}(),
                    Dict{String,Int}(),
                    String[],
                    A.nrdv,
                    A.npdv,
                    A.npv,
                    A.nlv,
                    ev_pol_locs,
                    ev_diff_pols_loc,
                    A.nomul,
                    order(A),
                    new_ctx,
                    A.inp)
end

function evaluate_parameter(pol :: OrePoly, point :: Int, A :: OreAlg) 
    T = eltype_co(A)
    ncoeffs = [T(evaluate(c,[point]).data) for c in coeffs(pol)]
    return OrePoly(ncoeffs,deepcopy(mons(pol)))
end
    
function evaluate_parameter(g :: Vector{OrePoly{T,M}}, point :: Int, A :: OreAlg)  where {T, M}
    return [evaluate_parameter(a, point, A) for a in g]
end
    
    
function evaluate_parameter(g :: Vector{Vector{OrePoly{T,M}}}, point :: Int,  A :: OreAlg) where  {T,M}
    return [evaluate_parameter(a, point, A) for a in g]
end
    
function evaluate_parameter(a :: Any,:: Int,:: OreAlg)
    return a 
end

function evaluate_parameter(p ::Int, A :: OreAlg, args...)
    return (evaluate_parameter(arg,p,A) for arg in args)
end

