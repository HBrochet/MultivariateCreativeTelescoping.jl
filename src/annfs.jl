
# function annfs_init(p :: OrePoly, A :: OreAlg)
#     inp = A.inp

#     pdv1,pdv2 = inp.poldiffvars
#     push!(pdv1,"_t")
#     push!(pdv2,"_dt") 
#     inp.poldiffvars = (pdv1,pdv2)

#     push!(inp.polvars,"_u")
#     push!(inp.polvars,"_v")
#     ord = 
#     inp.order = "lex u v "
# end
function annfs(p :: OrePoly, A :: OreAlg)
    gens = annfs_gens(p,A)
    gb = f5(gens,A)
    gb = annfs_rem_uv!(gb,A)
    return annfs_subsitute(gb,A)
end

function annfs_gens(p :: OrePoly{T,M},A :: OreAlg) where {T,M}
    t = parse_OrePoly("t",A)
    u = parse_OrePoly("u",A)
    v = parse_OrePoly("v",A)
    gens = OrePoly{T,M}[]
    
    push!(gens,parse_OrePoly("u*v-1",A))
    tmp = sub(p,parse_OrePoly("t*u",A),A)
    push!(gens,tmp)
    for i in 1:A.npdv-1
        tmp = mul(diff(p,i,A),parse_OrePoly("v*dt",A),A)
        tmp = add!(tmp,OrePoly([one(ctx(A))],[makemon(i,A)]),A)
        push!(gens,tmp)
    end
    return gens 
end

function annfs_rem_uv!(gb :: Vector{OrePoly{T,M}},A :: OreAlg) where {T,M} 
    t = Int[] 
    N = nvars(A)
    for i in 1:length(gb)
        m = lm(gb[i])
        if m[N-1] > 0 || m[N-2] > 0 
            push!(t,i)
        end
    end
    deleteat!(gb,t)
    return gb 
end


function annfs_subsitute(gb :: Vector{OrePoly{T,M}},A :: OreAlg) where{T,M} 
    tdt = makemon(A.npdv,A)*makemon(2*A.npdv,A)
    res = OrePoly{T,M}[]
    for g in gb
        degt = maximum(mon(g,i)[2*A.npdv] for i in 1:length(g))
        degdt = maximum(mon(g,i)[A.npdv] for i in 1:length(g))
        if degt > degdt 
            g = mul(parse_OrePoly("dt^($(degt-degdt))",A),g,A)
        elseif degdt > degt 
            g = mul(parse_OrePoly("t^($(degdt-degt))",A),g,A)
        end
        sub = parse_OrePoly("-s-1",A)
        tmp = zero(A)
        for (c,m) in g 
            @assert m[A.npdv] == m[A.npdv*2]
            d =  m[A.npdv]
            p = OrePoly([c],[m/tdt^d])
            for i in 1:d 
                p = mul(p,sub,A)
            end
            tmp = add!(tmp,p,A)
        end
        normalize!(tmp,A)
        push!(res,tmp)
    end
    return res 
end

function annfs_fglm(gb2 :: Vector{OrePoly{T,M}},A :: OreAlg) where{T,M} 
    s = parse_OrePoly("s",A)
    basis = OrePoly{T,M}[]
    p = one(A)
    while true
     # à continuer quand j'aurai refais une passe sur MCT 
    end 
end

