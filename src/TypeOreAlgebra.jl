abstract type AbsOreAlgebra{C,M,O} end;

# used to easily change from algebra 
mutable struct OreAlgInput
    order :: String 
    char :: Int 
    ratvars :: Vector{String}
    ratdiffvars :: Tuple{Vector{String},Vector{String}}
    poldiffvars :: Tuple{Vector{String},Vector{String}}
    polvars :: Vector{String}
    locvars :: Tuple{Vector{String},Vector{String}}
    nomul :: Vector{String}
end

struct OreAlg{T, Tbuf, C <:AbsContextCoeff{T,Tbuf}, M <: AbsOreMonomial, O <: AbsMonomialOrder} <: AbsOreAlgebra{C,M,O}
    strvar_to_indexp :: Dict{String,Int} # map strings of vars represented in monomials to corresponding index in the exponent
    indexp_to_strvar :: Vector{String} # same but in the other direction
    ratvars :: Dict{String,T}  # maps strings or ratvars represented in the coeff to an object of type T 
    rvars_to_int :: Dict{String,Int} # maps strings or ratvars represented in the coeff to their associated diff var if they have one 
    int_to_rvars :: Vector{String} # same but in the other direction


    nrdv :: Int # nb of pair (rat var, diff var) 
    npdv :: Int  # nb of pair of (pol var, diff var)
    npv :: Int # nb of pol vars without associated diff var
    nlv :: Int # nb of loc vars
    # order of the variables: rdv, diff pdv, pol pdv, pv, lv
    
    pols_loc :: Vector{OrePoly{T,M}} # pol associated to the loc var
    diff_pols_loc :: Vector{Vector{OrePoly{T,M}}} #derivatives of pols_loc
    nomul :: Vector{Int} # indices of the variable with which we can't multiply
    ctx :: C 
    inp :: OreAlgInput
end

"""
    OreAlg(;char :: Int = 0,
    ratvars :: Vector{String} = String[],
    ratdiffvars :: Tuple{Vector{String},Vector{String}}=(String[],String[]), 
    poldiffvars :: Tuple{Vector{String},Vector{String}}=(String[],String[]), 
    polvars :: Vector{String} = String[],
    locvars :: Tuple{Vector{String},Vector{String}} = (String[],String[]),
    order :: String = "",
    nomul :: Vector{String} = String[]

For more information see the online documentation
"""
function OreAlg(;order :: String = "",
    char :: Int = 0,
    ratvars :: Vector{String} = String[],
    ratdiffvars :: Tuple{Vector{String},Vector{String}}=(String[],String[]), 
    poldiffvars :: Tuple{Vector{String},Vector{String}}=(String[],String[]), 
    polvars :: Vector{String} = String[],
    locvars :: Tuple{Vector{String},Vector{String}} = (String[],String[]),
    nomul :: Vector{String} = String[]
        )
    ### Create context for coefficients
    ctx = make_coeff_context(char,ratvars,ratdiffvars)

    ### Initialising variable name tables
    sti, its, rv, rti, itr = make_var_tables(ratvars,ratdiffvars,poldiffvars,polvars,locvars,ctx)

    # Type of the monomials 
    M = OreMonVE{length(its),Int16}
    #Type of coefficients
    T = eltype1_ctx(ctx)

    ### Creating the monomial order
    ord = make_order(order,sti,Val(M))

    inp = OreAlgInput(order, char,ratvars, ratdiffvars,poldiffvars,polvars,locvars,nomul)

    ### Computing the polynomials w.r.t. which we localize and their derivatives
    tmpA = OreAlg{eltype1_ctx(ctx),eltype2_ctx(ctx), typeof(ctx),M,typeof(ord)}(sti,
                                            its,
                                            rv,
                                            rti,
                                            itr,
                                            length(ratdiffvars[2]),
                                            length(poldiffvars[1]),
                                            length(polvars),
                                            length(locvars[1]),
                                            OrePoly{T,M}[],
                                            Vector{OrePoly{T,M}}[],
                                            Int[],
                                            ctx,
                                            inp)
    
    (polsloc, diff_pols_loc) = compute_locpol_and_derivs(locvars[2],tmpA)

    nomul_ = [sti[s] for s in nomul]
    return OreAlg{eltype1_ctx(ctx),eltype2_ctx(ctx), typeof(ctx),M,typeof(ord)}(sti,
                                            its,
                                            rv,
                                            rti,
                                            itr,
                                            length(ratdiffvars[2]),
                                            length(poldiffvars[1]),
                                            length(polvars),
                                            length(locvars[1]),
                                            polsloc,
                                            diff_pols_loc,
                                            nomul_,
                                            ctx,
                                            inp)
end


function make_coeff_context(char :: Int, ratvars :: Vector{String}, ratdiffvars :: Tuple{Vector{String},Vector{String}})
    # coefficients will be rational functions 
    if (length(ratvars) + length(ratdiffvars[1])) > 0 
        if char > 0 
            S,vars = polynomial_ring(Native.GF(char),vcat(ratdiffvars[1],ratvars))
            F = fraction_field(S)
            return RatFunCtx(F,S,vars,char)
        else 
            S,vars = polynomial_ring(ZZ,vcat(ratdiffvars[1],ratvars))
            F = fraction_field(S)
            return RatFunQQCtx(F,S,vars)
        end
    # coefficients with be rationals or int mod p 
    else 
        if char > 0 
            return Nmod32Γ(char)
        else 
            return QQCtx()
        end
    end
end

function make_var_tables(ratvars :: Vector{String},
    ratdiffvars :: Tuple{Vector{String},Vector{String}}, 
    poldiffvars :: Tuple{Vector{String},Vector{String}}, 
    polvars :: Vector{String},
    locvars :: Tuple{Vector{String},Vector{String}},
    ctx :: AbsContextCoeff)


    T = eltype1_ctx(ctx)
    sti = Dict{String,Int}()
    its = String[]
    rv = Dict{String,T}()
    rti = Dict{String,Int}()
    itr = String[]

    for (i,s) in enumerate(ratvars)
        rv[s] = ctx.F(ctx.vars[i+length(ratdiffvars[1])])
    end

    for (i,s) in enumerate(ratdiffvars[2])
        sti[s] = i
        push!(its,s)
    end
    for (i,s) in enumerate(ratdiffvars[1])
        rv[s] = ctx.F(ctx.vars[i])
        push!(itr,s)
        rti[s] = i
    end

    it = length(ratdiffvars[2])
    for l in 2:-1:1 
        for i in 1:length(poldiffvars[1])
            sti[poldiffvars[l][i]] = it + i
            push!(its,poldiffvars[l][i])
        end
        it = it + length(poldiffvars[1])
    end
    for (i,s) in enumerate(polvars)
        sti[s] = it + i
        push!(its,s)
    end
    it = it + length(polvars)

    for (i,s) in enumerate(locvars[1])
        sti[s] = it + i 
        push!(its,s)
    end
    return sti, its, rv, rti, itr
end

### Initialize polsloc and diff_pols_loc
function compute_locpol_and_derivs(pollocvars ::Vector{String},A ::OreAlg)
    T = eltype1_ctx(ctx(A))
    M = eltype_mo(A)
    polsloc = OrePoly{T,M}[]
    diff_pols_loc = Vector{OrePoly{T,M}}[]
    for (i,str) in enumerate(pollocvars)
        push!(diff_pols_loc, OrePoly{T,M}[])
        pol = parse_OrePoly(str,A)
        push!(polsloc,pol)
        append!(diff_pols_loc[i],[diff(pol,i,A) for i in 1:A.nrdv])
        append!(diff_pols_loc[i],[diff(pol,i + A.nrdv+A.npdv,A) for i in 1:A.npdv])
    end
    return (polsloc,diff_pols_loc) 
end

Base.show(io::IO, A::OreAlg) = print(io,"OreAlgebra")





function OreAlg(inp :: OreAlgInput) 
    return OreAlg(order = inp.order;
                char = inp.char,
                ratvars = inp.ratvars,
                ratdiffvars = inp.ratdiffvars, 
                poldiffvars = inp.poldiffvars, 
                polvars = inp.polvars,
                locvars = inp.locvars,
                nomul = inp.nomul)
end



order(A :: OreAlg{T,Tbuf,K,M,O}) where {T,Tbuf,K,M,O} = O()
ctx(A :: OreAlg) = A.ctx
eltype_co(A :: OreAlg{T,Tbuf,C,M,O}) where {T, Tbuf, C <: AbsContextCoeff{T, Tbuf},M,O} = T
eltype_mo(A :: OreAlg{T,Tbuf,K,M,O}) where {T,Tbuf,K,M,O} = M
makemon(i :: Integer, A :: OreAlg{T, TT, K, M, O}) where {T, TT, K,O,N,E,M <: OreMonVE{N,E}} = OreMonVE(SVector{N,E}(i==j ? E(1) : E(0) for j in 1:N))

function divide(m1 :: OreMonVE, m2 :: OreMonVE{N,E}, A :: AbsOreAlgebra) where {N,E}
    for i in A.nomul
        if m1[i] != m2[i]
            return false
        end
    end
    for i in 1:N 
        if m1[i] > m2[i]
            return false 
        end
    end
    return true 
end

function iscompatible(m1 :: OreMonVE, m2 :: OreMonVE, A :: AbsOreAlgebra) 
    for i in A.nomul 
        if m1[i] != m2[i] 
            return false 
        end
    end 
    return true
end

function maxdeg(p :: OrePoly,j :: Int)
    return maximum(mon(p,i)[j] for i in 1:length(p))
end  

function maxdeg(g :: Vector{OrePoly{K,M}},j :: Int) where {K, M}
    return maximum(maxdeg(p,j) for p in g)
end  


makepoly(c :: K, m :: M) where {K, M <: AbsOreMonomial} = OrePoly(K[c],M[m]) 
undefOrePoly(n :: Integer,A :: OreAlg{T,Tbuf,K,M,O}) where {T,Tbuf,K<:AbsContextCoeff{T,Tbuf},M,O} = OrePoly(Vector{T}(undef,n),Vector{M}(undef,n))

Base.zero(A :: OreAlg{T, Tbuf, K,M,O}) where {T,Tbuf,K<:AbsContextCoeff{T,Tbuf},M,O} = OrePoly(T[],M[])
Base.one(A :: OreAlg{T,Tbuf,K,M,O}) where {T,Tbuf,K<:AbsContextCoeff{T,Tbuf},M,O} = OrePoly(T[convert(1,ctx(A))],M[makemon(-1,A)])
nvars(A :: OreAlg{T,Tbuf,K,M,O}) where {T,Tbuf,N,E,K,M <:OreMonVE{N,E},O} = N


function diff(pol_ :: OrePoly, ind_ :: Integer, A ::AbsOreAlgebra)
    pol = deepcopy(pol_)
    ind = ind_
    if ind <= A.nrdv 
        cos = coeffs(pol)
        for i in 1:length(cos)
            cos[i] = derivative(cos[i], ctx(A).vars[ind])
        end
        normalize!(pol,A)
        return pol 
    end
    if ind <= A.nrdv + A.npdv 
        ind += A.npdv 
    end
    for j in 1:length(pol)
        if mon(pol,j)[ind] > 0
            co = mul(coeffs(pol)[j],convertn(mon(pol,j)[ind],ctx(A)),ctx(A))
            mo = mons(pol)[j] / makemon(ind,A)
            pol[j] = (co,mo)
        else 
            coeffs(pol)[j] = zero(ctx(A))
        end 

    end
    normalize!(pol,A)
    return pol
end


function denominator(p :: OrePoly{T,M}, A :: OreAlg) where {M, T <: Generic.FracFieldElem}
    if length(p) == 0 
        return ctx(A).R(1)
    end
    res = ctx(A).R(1)
    for c in coeffs(p)
        res = lcm(res, Nemo.denominator(c))
    end
    return res 
end



function isholonomic(gb :: Vector{OrePoly{T,M}},A :: Alg) where {T,M, Alg <:AbsOreAlgebra}
    #println("testing holonomy")
    N = nvars(A)
    set = [i for i in 1:A.nrdv+2*A.npdv]
    n = A.npdv
    # returns an iterator over every subset of set of size n+1
    it = powerset(set,n+1,n+1) 

    for subset in it 
        gfound = false
        for i in 1:length(gb)
            if A.nlv > 0 && mon(gb[i],1)[A.nrdv+2*A.npdv+A.npv+1] > 0
                continue
            end 
            # test if gb[i] verifies the criterion
            b = true  
            for j in 1:N
                if !(j in subset) && mon(gb[i],1)[j] != 0
                    b = false
                    break
                end
            end
            if b 
                gfound = true
                break
            end
        end
        if !gfound 
            #println(subset)
            return false
        end
    end
    return true
end

# check if p and q commute
function commute(p :: OrePoly, q :: OrePoly, A :: OreAlg)
    if A.nlv > 0 
        return false
    end

    # compute the support of p and q (variable i appear in p iff ep[i] > 0)
    ep = mon(p,1)
    for i in 2:length(p)
        ep = lcm(ep,mon(p,i))
    end
    eq = mon(q,1)
    for i in 2:length(q)
        eq = lcm(eq,mon(q,i))
    end

    for i in A.nrdv+1:A.nrdv+A.npdv
        if (ep[i] > 0 && eq[i+A.npdv] > 0) || (ep[i+A.npdv] > 0 && eq[i] > 0 )
            return false 
        end
    end

    if A.nrdv == 0
        return true
    end
    for i in 1:A.nrdv 
        v = numerator(A.ratvars[A.int_to_rvars[i]])
        if eq[i] > 0 
            for j in 1:length(p) 
                if v in vars(coeff(p,j))
                    return false
                end
            end
        end
        if ep[i] > 0 
            for j in 1:length(q)
                if v in vars(coeff(q,j))
                    return false 
                end
            end
        end
    end
    return true
end

function leadingterms(v :: Vector{OrePoly{T,M}}) where {T,M}
    return [OrePoly([coeff(p,1)],[mon(p,1)]) for p in v]
end