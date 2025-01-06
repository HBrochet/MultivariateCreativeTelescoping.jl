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
    if ctx(A) isa UnivRatFunCtx
        nctx = UnivRatFunModpCtx(F,R,nvars,prime)
        ev_ratvars = Dict{String,UnivRatFunModp}()
        T = UnivRatFunModp
    else
        nctx = RatFunModpCtx(F,R,nvars,prime)
        ev_ratvars = Dict{String,RatFunModp}()
        T = RatFunModp
    end
    for k in keys(A.ratvars)
        ev_ratvars[k] = F(change_coefficient_ring(CRing, numerator(A.ratvars[k], false)))
    end
    inp = deepcopy(A.inp)
    inp.char = prime

    M = eltype_mo(A)
    tmpA = OreAlg{eltype1_ctx(nctx),typeof(nctx),eltype_mo(A),eltype_ord(A)}(A.strvar_to_indexp,
                                                                A.indexp_to_strvar,
                                                                ev_ratvars,                                                               
                                                                A.drvars_to_int,
                                                                A.int_to_drvars,
                                                                A.nrdv,
                                                                A.npdv,
                                                                A.npv,
                                                                A.nlv,
                                                                OrePoly{T,M}[],
                                                                Vector{OrePoly{T,M}}[],
                                                                A.nomul,
                                                                A.order,
                                                                nctx,
                                                                A.varord,
                                                                A.inp)

    pl = [change_coefficient_field(CRing,p,tmpA) for p in A.pols_loc]
    dpl = [change_coefficient_field(CRing,v,tmpA) for v in A.diff_pols_loc]

    return OreAlg{eltype1_ctx(nctx),typeof(nctx),eltype_mo(A),eltype_ord(A)}(A.strvar_to_indexp,
                                                                A.indexp_to_strvar,
                                                                ev_ratvars,                                                               
                                                                A.drvars_to_int,
                                                                A.int_to_drvars,
                                                                A.nrdv,
                                                                A.npdv,
                                                                A.npv,
                                                                A.nlv,
                                                                pl,
                                                                dpl,
                                                                A.nomul,
                                                                A.order,
                                                                nctx,
                                                                A.varord,
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
    tmpA = OreAlg{eltype1_ctx(nctx),typeof(nctx),eltype_mo(A),eltype_ord(A)}(A.strvar_to_indexp,
                                                                A.indexp_to_strvar,
                                                                nratvars,                                                               
                                                                A.drvars_to_int,
                                                                A.int_to_drvars,
                                                                A.nrdv,
                                                                A.npdv,
                                                                A.npv,
                                                                A.nlv,
                                                                OrePoly{T,M}[],
                                                                Vector{OrePoly{T,M}}[],
                                                                A.nomul,
                                                                A.order,
                                                                nctx,
                                                                A.varord,
                                                                A.inp)

    pl = [change_coefficient_field(p,tmpA) for p in A.pols_loc]
    dpl = [change_coefficient_field(v,tmpA) for v in A.diff_pols_loc]

    return OreAlg{eltype1_ctx(nctx),typeof(nctx),eltype_mo(A),eltype_ord(A)}(A.strvar_to_indexp,
                                                                A.indexp_to_strvar,
                                                                nratvars,                                                               
                                                                A.drvars_to_int,
                                                                A.int_to_drvars,
                                                                A.nrdv,
                                                                A.npdv,
                                                                A.npv,
                                                                A.nlv,
                                                                pl,
                                                                dpl,
                                                                A.nomul,
                                                                A.order,
                                                                nctx,
                                                                A.varord,
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
    num = numerator(p, false)
    den = Nemo.denominator(p, false)
    F = ctx(A).F
    return F(change_coefficient_ring(R,num)) / F(change_coefficient_ring(R,den))
end

function change_coefficient_field(R :: Ring, p :: Generic.FracFieldElem{ZZPolyRingElem}, A :: OreAlg)
    num = numerator(p, false)
    den = Nemo.denominator(p, false)
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
    tmpA = OreAlg{T,typeof(new_ctx),eltype_mo(A),eltype_ord(A)}(A.strvar_to_indexp,
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
                    A.varord,
                    A.inp)
    ev_pol_locs = evaluate_parameter(A.pols_loc, p, tmpA)
    ev_diff_pols_loc = evaluate_parameter(A.diff_pols_loc, p, tmpA)
    return OreAlg{eltype1_ctx(new_ctx),typeof(new_ctx),eltype_mo(A),eltype_ord(A)}(A.strvar_to_indexp,
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
                    A.varord,
                    A.inp)
end

function evaluate_parameter(pol :: OrePoly, point :: Int, A :: OreAlg) 
    T = eltype_co(A)
    cos = coeffs(pol)
    len = length(cos)
    ncoeffs = Vector{T}(undef,len)
    for i in 1:len
        ncoeffs[i] =convertn(evaluate(cos[i],point).data,ctx(A))
    end
    return OrePoly(ncoeffs,deepcopy(mons(pol)))
end
    
function evaluate_parameter(g :: Vector{OrePoly{T,M}}, point :: Int, A :: OreAlg)  where {T, M}
    return [evaluate_parameter(a, point, A) for a in g]
end

function evaluate_parameter(d :: Dict{M,OrePoly{T,M}}, point :: Int, A :: OreAlg)  where {T, M}
    return Dict(k=>evaluate_parameter(d[k],point,A) for k in keys(d))
end
    
    
function evaluate_parameter(g :: Vector{Vector{OrePoly{T,M}}}, point :: Int,  A :: OreAlg) where  {T,M}
    return [evaluate_parameter(a, point, A) for a in g]
end
   
function evaluate_parameter(t :: Tuple,p:: Int,A:: OreAlg)
    return (evaluate_parameter(t[1],p,A), evaluate_parameter(t[2],p,A))
end
function evaluate_parameter(a :: Any,:: Int,:: OreAlg)
    return a 
end

function evaluate_parameter(p ::Int, A :: OreAlg, args...)
    return (evaluate_parameter(arg,p,A) for arg in args)
end

function evaluate_parameter_cbl(cbl :: UnivRatFunModp ,randpoints :: Vector{Int},A :: OreAlg)
    return [UInt32(Nemo.evaluate(Nemo.denominator(cbl, false),point).data) for point in randpoints]
end

### Evaluate the single parameter at many points using multievaluation 

# Type unstable !!
function evaluate_parameter_many(glen :: Int, randpoints :: Vector{Int}, nA :: OreAlg, args...;denisone ::Val{TT}= Val(false)) where TT
    vpoints = new_rand_points(randpoints,Int(char(nA)),glen)
    _vpoints = [UInt(p) for p in vpoints]
    t =  ntuple(length(args)) do i
        evaluate_parameter_many(args[i],_vpoints,nA;denisone=Val(TT))
    end
    return vpoints, t
end


function evaluate_parameter_many(glen :: Int, randpoints :: Vector{Int}, nA :: OreAlg, a,b;denisone ::Val{T}= Val(false)) where T
    vpoints = new_rand_points(randpoints,Int(char(nA)),glen)
    _vpoints = [UInt(p) for p in vpoints]
    a = evaluate_parameter_many(a,_vpoints,nA;denisone=Val(T))
    b = evaluate_parameter_many(b,_vpoints,nA;denisone=Val(T))

    return vpoints, tuple(a,b)
end

function evaluate_parameter_many(glen :: Int, randpoints :: Vector{Int}, nA :: OreAlg, a,b,c;denisone ::Val{T}= Val(false)) where T
    vpoints = new_rand_points(randpoints,Int(char(nA)),glen)
    _vpoints = [UInt(p) for p in vpoints]
    a = evaluate_parameter_many(a,_vpoints,nA;denisone=Val(T))
    b = evaluate_parameter_many(b,_vpoints,nA;denisone=Val(T))
    c = evaluate_parameter_many(c,_vpoints,nA;denisone=Val(T))
    return vpoints, tuple(a,b,c)
end


function new_rand_points(vec :: Vector{Int}, p :: Int,len :: Int)
    l = 0 
    res = Int[]
    while l<len  
        point = mod(rand(Int),p)
        if !(point in vec) && !(point in res)
            push!(res,point)
            l += 1
        end
    end
    return res
end


function evaluate_parameter_many(pol :: OrePoly, v :: Vector{UInt}, nA:: OreAlg; denisone :: Val{B} = Val(false)) where B
    T = eltype_co(nA)
    M = eltype_mo(nA)
    if B 
        tmp = Vector{Vector{UInt}}(undef,length(pol))
    else
        tmp = Vector{Tuple{Vector{UInt}, Vector{UInt}}}(undef,length(pol))
    end
    for i in 1:length(pol)
        tmp[i] = evaluate_many(coeff(pol,i),v,denisone = Val(B))
    end
    res = Vector{OrePoly{T,M}}(undef,length(v))
    for l in 1:length(v)
        cs = Vector{T}(undef,length(pol))
        for i in 1:length(pol)
            if B 
                cs[i] = convertn(tmp[i][1][l],ctx(nA))
            else
                cs[i] = mul(T(tmp[i][1][l]),inv(T(tmp[i][2][l]),ctx(nA)),ctx(nA))
            end
        end
        res[l] = OrePoly(cs,deepcopy(mons(pol)))
    end
    return res
end


    
function evaluate_parameter_many(g :: Vector{OrePoly{TT,M}}, v :: Vector{UInt}, nA :: OreAlg;denisone ::Val{B} = Val(false))  where {B, M,TT}
    T = eltype_co(nA)
    tmp = Vector{Vector{OrePoly{T,M}}}(undef,length(g))
    for i in 1:length(g)
        tmp[i] = evaluate_parameter_many(g[i],v,nA,denisone=Val(B))
    end
    res = Vector{Vector{OrePoly{T,M}}}(undef,length(v))
    for l in 1:length(v)
        gg = Vector{OrePoly{T,M}}(undef,length(g))
        for i in 1:length(g)
            gg[i] = tmp[i][l]
        end
        res[l] = gg
    end
    return res
end

function evaluate_parameter_many(d :: Dict{M,OrePoly{TT,M}}, v :: Vector{UInt}, nA :: OreAlg;denisone ::Val{B}= Val(false))  where {TT, M,B}
    T = eltype_co(nA)
    tmp = Dict{M,Vector{OrePoly{T,M}}}()
    for k in keys(d)
        tmp[k] = evaluate_parameter_many(d[k],v,nA,denisone=Val(B))
    end
    res = Vector{Dict{M,OrePoly{T,M}}}(undef,length(v))
    for l in 1:length(v)
        dd = Dict{M,OrePoly{T,M}}()
        for k in keys(d)
            dd[k] = tmp[k][l]
        end
        res[l] = dd
    end
    return res
end

#rendu ici 
function evaluate_parameter_many(mat :: Generic.MatSpaceElem{Generic.FracFieldElem{Nemo.fpPolyRingElem}}, v :: Vector{UInt}, nA :: OreAlg;denisone ::Val{T}= Val(false)) where T
    nc = number_of_columns(mat)
    nr = number_of_rows(mat)

    if T 
        tmp = Matrix{Vector{UInt}}(undef,nr,nc)
    else
        tmp = Matrix{Tuple{Vector{UInt}, Vector{UInt}}}(undef,nr,nc)
    end
    for j in 1:nc
        for i in 1:nr 
            tmp[i,j] = evaluate_many(mat[i,j],v;denisone=Val(T))
        end
    end

    S = Native.GF(Int(char(nA)))
    M = matrix_space(S,nr,nc)
    res = Vector{fpMatrix}(undef,length(v))

    for l in 1:length(v)
        evmat = M()
        for j in 1:nc
            for i in 1:nr 
                if T 
                    evmat[i,j] = S(tmp[i,j][l])
                else
                    evmat[i,j] = S(tmp[i,j][1][l]) //  S(tmp[i,j][2][l])
                end 
            end
        end
        res[l] = evmat
    end
    return res
end
function evaluate_parameter_many(p :: MCTParam, v :: Vector{UInt}, nA :: OreAlg;denisone ::Val{T}= Val(false)) where T
    return [p for i in 1:length(v)] 
end

# function evaluate_parameter(g :: Vector{Vector{OrePoly{T,M}}}, point :: Int,  A :: OreAlg) where  {T,M}
#     return [evaluate_parameter(a, point, A) for a in g]
# end
   
# function evaluate_parameter(t :: Tuple,p:: Int,A:: OreAlg)
#     return (evaluate_parameter(t[1],p,A), evaluate_parameter(t[2],p,A))
# end
# function evaluate_parameter(a :: Any,:: Int,:: OreAlg)
#     return a 
# end

