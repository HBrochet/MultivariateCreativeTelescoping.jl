struct WCParam{T, S} 
    gb_param :: T 
    bound :: Int
    morestep :: Int
end 

function wc_param(;method ::Val{F} = Val(:f4),
    bound :: Int = 1,
    morestep :: Int  = 0,
    geobucket ::Val{B} = Val(false),
    stophol  :: Val{C} = Val(false),
    stat :: Val{D} = Val(false),
    debug :: Val{E} = Val(false),
    select_reducer :: Val{G} = Val(:shortest),
    tracer :: Val{H} = Val(:none)
    ) where {B,C,D,E,F,G,H} 
    if F == :f5 
        param =  F5Param{B,C,D,E}()
    elseif (F == :f4) || (F == :buchberger) 
        param = F4Param{B,C,D,E,G,H}()
    end
    return WCParam{typeof(param),F}(param,bound,morestep) 
end

method(:: WCParam{T,S}) where {T,S} = S
stophol(p:: WCParam) = stophol(p.gb_param)
stat(p:: WCParam) = stat(p.gb_param)
debug(p:: WCParam)  = debug(p.gb_param)
# the two functions below check whether a tracer is used
learn(p:: WCParam)  = learn(p.gb_param)
apply(p:: WCParam)  = apply(p.gb_param)

function wc_param_learn_to_apply(p:: WCParam) 
    param = f4param_learn_to_apply(p.gb_param)
    return WCParam{typeof(param), :f4}(param,p.bound,p.morestep)
end


struct WCTrace 
    v :: Vector{F4Trace}  
end
wc_trace() = WCTrace(F4Trace[])




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
function saturation(g :: Vector{OrePoly{T,M}}, i :: Int,  A :: OreAlg, trace ::WCTrace ; param :: WCParam = wc_param()) where {T,M}
    bound = param.bound # bound on the degree in the loc variable
    morestep = param.morestep
    ctr = 1 
    varT = OrePoly([one(ctx(A))],[makemon(nvars(A),A)])
    tmp = mul(A.pols_loc[1],varT,A)
    tmp = sub(tmp,one(A),A)
    push!(g, tmp)
    

    if method(param) == :buchberger
        gb = Buchberger2(g,A,param = param.gb_param)
    elseif method(param) == :f4 
        if learn(param)
            gb ,tr = f4(g,A,param = param.gb_param)
            push!(trace.v, tr)
        elseif apply(param)
            gb = f4(g,A,param = param.gb_param,trace = trace.v[ctr])
        else
            gb = f4(g,A,param = param.gb_param)
        end
    else
        gb = f5(g,A,param = param.gb_param)
    end

    gb, lastpow = init_gb_sat(gb, bound,i,A)
    # println("gb")
    # prettyprint(gb,A)
    res = filter(p -> mon(p,1)[i] == 0, gb) 
    debug(param) && @debug "gb initialized, it contains $(length(gb)) vectors"
    hol = isholonomic(res,A) 
    while !hol || morestep > 0 
        bound = bound + 1
        if hol 
            morestep -= 1 
        end
        update_gb_sat!(gb,lastpow,bound,i,A)
        ctr += 1 
        # ctr == 5 && error("fin")
        if method(param) == :buchberger
            gb = Buchberger2(gb,A,param = param.gb_param)
        elseif method(param) == :f4
            if learn(param)
                gb ,tr = f4(gb,A,param = param.gb_param)
                push!(trace.v, tr)
            elseif apply(param)
                gb = f4(gb,A,param = param.gb_param,trace = trace.v[ctr])
            else
                gb = f4(gb,A,param = param.gb_param)
            end
        else
            gb = f5(gb,A,param = param.gb_param)
        end
        res = filter(p -> mon(p,1)[i] == 0, gb) 
        hol = isholonomic(res,A)
    end
    if learn(param)
        return res, trace
    else
        return res 
    end
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


function weyl_closure_internal(gens::Vector{OrePoly{T,M}}, A :: OreAlg, init :: WeylClosureInit, param :: WCParam; trace :: WCTrace = wc_trace()) where {T,M}
    nA = init.A
    sl = singular_locus(gens,A,init.SLI)
    if sl == one(A)
        return f4(gens,A)
    else 
        nemo_p = prod(fact[1] for fact in factor(sl) )
        p = OrePoly([nemo_p],[makemon(-1,A)])
    end
    # change algebra so that _T corresponds to 1/p
    (polsloc,diff_pols_loc) = compute_locpol_and_derivs([mystring(p,A)],nA)
    empty!(nA.pols_loc)
    append!(nA.pols_loc,polsloc)
    empty!(nA.diff_pols_loc) 
    append!(nA.diff_pols_loc, diff_pols_loc) 
    ngens = map_algebras(gens,A,nA)
    if learn(param)
        gb, tr = saturation(ngens,nA.nrdv+2*nA.npdv +nA.npv+nA.nlv,nA,trace, param = param) 
        return map_algebras(gb,nA,A), tr 
    else
        gb = saturation(ngens,nA.nrdv+2*nA.npdv +nA.npv+nA.nlv,nA,trace, param = param) 
        return map_algebras(gb,nA,A)
    end
end



# wrapper for using crt with tracer = Val(:apply) 
function weyl_closure_internal_crt(trace :: WCTrace, gens::Vector{OrePoly{T,M}}, A :: OreAlg, init :: WeylClosureInit, param :: WCParam) where {T,M}
    par = wc_param_learn_to_apply(param)
    # println("learn $(learn(par)) et $(apply(par))")
    return weyl_closure_internal(gens, A, init, par, trace = trace)
end

# wrapper for using crt with tracer = Val(:learn) 
function weyl_closure_internal_crt(gens::Vector{OrePoly{T,M}}, A :: OreAlg, init :: WeylClosureInit, param :: WCParam; tracer :: Val{B} = Val(false)) where {T,M,B}
    # println("learn $(learn(param)) et $(apply(param))")
    B && @assert learn(param) # for compatibility with compute_with_crt
    return weyl_closure_internal(gens, A, init, param, trace = wc_trace())

end



"""
    weyl_closure(gens::Vector{OrePoly{T,M}}, A :: OreAlg, init :: WeylClosureInit)

Return a Gröbner basis of a submodule of the Weyl closure of the module generated by gens.
"""
function weyl_closure(gens::Vector{OrePoly{T,M}}, A :: OreAlg, init :: WeylClosureInit; param :: WCParam = wc_param()) where {T,M}
    if ctx(A) isa UnivRatFunModpCtx || ctx(A) isa RatFunModpCtx || ctx(A) isa UnivRatFunModPCtx || ctx(A) isa RatFunModPCtx || ctx(A) isa NmodΓ
        return weyl_closure_internal(gens,A, init,param)
    else 
        par = crt_param(comp=Val(:slow),tracer = Val(learn(param)))
        return compute_with_CRT(weyl_closure_internal_crt,A, gens, A, init,param;param = par)
    end
end



function deg_d(m :: OreMonVE, A :: OreAlg) 
    return sum(m[i] for i in 1:A.nrdv + A.npdv)
end

function initial_form(gens :: Vector{OrePoly{T,K}}, A :: OreAlg) where {T,K}
    res = OrePoly{T,K}[]
    N = nvars(A) 
    E = Int16
    for g in gens 
        c = deg_d(mon(g,1),A) 
        l= 1
        while l+1 <= length(g)
            deg_d(mon(g,l+1),A) < c && break
            l +=1 
        end
        p = undefOrePoly(l,A)
        for i in 1:l 
            p.coeffs[i] = coeff(g,i)
            m = mon(g,i)
            m2 = OreMonVE(SVector{N,E}(i <= A.nrdv + A.npdv ? m[i] : 0 for i in 1:N))
            m3 = OreMonVE(SVector{N,E}(((i > A.nrdv + A.npdv + A.npv) && i<= A.nrdv + 2*A.npdv + A.npv) ? m[i-A.nrdv - A.npdv - A.npv] : 0 for i in 1:N))
            p.mons[i] = m/m2*m3
        end
        push!(res,p)
    end
    return res 
end

function intersect_with_cx!(gens :: Vector{OrePoly{T,K}},A ::OreAlg) where {T,K} 
    to_remove = Int[] 
    for (j,g) in enumerate(gens) 
        m = mon(g,1)
        if any(m[i] > 0 for i in A.nrdv + A.npdv + A.npv+1:A.nrdv + 2*A.npdv + A.npv)
            push!(to_remove,j)
        end
    end
    deleteat!(gens,to_remove)
    return gens 
end



function singular_locus2_init(gens :: Vector{OrePoly{T,K}},A ::OreAlg) where {T,K} 
    vars  = copy(A.inp.poldiffvars[2])
    append!(vars,A.inp.ratdiffvars[2])
    pvars = ["p$(v)" for v in vars]
    Tvars = ["T$(v)" for v in vars]

    order = "grevlex "*join(Tvars," ")*" > grevlex "*join(pvars," ")*" > grevlex " * join(vars," ") * " "*join(A.inp.ratdiffvars[2]," ")*" > grevlex "*join(A.inp.poldiffvars[1]," ") * " "*join(A.inp.polvars," ")
    inp = A.inp 
    inp.order = order 
    append!(inp.locvars[1],Tvars) 
    append!(inp.locvars[2],pvars)
    append!(inp.polvars,pvars)
    A2 = OreAlg(inp)
    nA = to_ore_alg_with_rat_coeffs(A2)
    return A2,nA 
end


function singular_locus2(gens :: Vector{OrePoly{T,K}},A ::OreAlg, A2 :: OreAlg, nA :: OreAlg) where {T,K} 
    vars  = copy(A.inp.poldiffvars[2])
    append!(vars,A.inp.ratdiffvars[2])

    ngens = parse_vector_OrePoly(mystring(gens,A),A2)
    gb = f4(ngens,A2)
    ini = initial_form(gb,A2)
    append!(ini, [parse_OrePoly("p$(v)*T$(v)-1",A2) for v in vars])
    gb2 = f4(ini,A2)
    intersect_with_cx!(gb2,A2)
    ngens = map_algebras(gb2,A2,nA)
    thelcm = lcm([Nemo.numerator(coeff(g,1)) for g in ngens])
    return thelcm
end


