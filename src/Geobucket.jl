const _d = 4 # the size of buckets in the geobucket grows as _d^i 

struct GeoBucket{K, M <: AbsOreMonomial} <: AbsOrePolynomial{K, M}
    tab1 :: Vector{OrePoly{K,M}}
    tab2 :: Vector{OrePoly{K,M}}
    active_tab_is_1 :: Vector{Bool} # contains the value true (resp. false) is tab1[i] is active (resp. tab2[i])
    indices :: Vector{Int} # index of the lm in the ith active tab 
end

# following Yan's paper we create an array tab s.t length(tab[i]) <= 4*_d^i 
# that store a polynomial of size at most _d^i
# It is stored in the first _d^i indices if indices[i] <= _d^i and between the indices 2*_d^i +1:end 
# otherwise. The remaining memory is preallocated for addition of polynomials 
# !!! Monomials are stored in increasing order in GeoBuckets

function GeoBucket(p :: OrePoly{K,M}) where {K, M<: AbsOreMonomial}
    lenp = length(p)
    l = int_log_upper(_d,lenp)
    tab1 = [OrePoly(Vector{K}(undef,2*_d^i),Vector{M}(undef,2*_d^i)) for i in 1:l]
    tab2 = [OrePoly(Vector{K}(undef,2*_d^i),Vector{M}(undef,2*_d^i)) for i in 1:l]
    active_tab_is_1 = [true for j in 1:l]
    indices = [0 for j in 1:l]
    for j in 1:lenp
        tab1[l][lenp-j+1] = p[j]
    end
    indices[l] = lenp
    return GeoBucket{K,M}(tab1,tab2,active_tab_is_1,indices) 
end

function init_GeoBucket!(g :: GeoBucket, p :: OrePoly)
    lenp = length(p)
    l = int_log_upper(_d,lenp)
    if l > length(g) 
        grow_to!(g,l)
    end
    tab = g.tab1[l]
    for j in 1:lenp
        tab[lenp-j+1] = p[j]
    end
    g.indices[l] = lenp
    g.active_tab_is_1[l] = true
    return g
end

tab1_is_active(g :: GeoBucket, i :: Int) = g.active_tab_is_1[i]
Base.length(g :: GeoBucket) = length(g.tab1)
Base.iszero(g :: GeoBucket, i ::Int) = g.indices[i] == 0
 
Base.show(io::IO, g :: GeoBucket) = print(io,"Geo bucket")
Base.show(io::IO, ::MIME"text/plain", g :: GeoBucket) = print(io,"Geo bucket")


function lm(g :: GeoBucket, i ::Int)
    if tab1_is_active(g,i)
        return mon(g.tab1[i],g.indices[i])
    else
        return mon(g.tab2[i],g.indices[i])
    end
end

function lc(g :: GeoBucket, i ::Int)
    if tab1_is_active(g,i)
        return coeff(g.tab1[i],g.indices[i])
    else
        return coeff(g.tab2[i],g.indices[i])
    end
end

function active_tab(g :: GeoBucket, i ::Int)
    if tab1_is_active(g,i)
        return g.tab1[i]
    else
        return g.tab2[i]
    end
end

function inactive_tab(g :: GeoBucket, i ::Int)
    if tab1_is_active(g,i)
        return g.tab2[i]
    else
        return g.tab1[i]
    end
end

function change_active_tab!(g :: GeoBucket,j :: Int)
    g.active_tab_is_1[j] = !g.active_tab_is_1[j] 
    nothing 
end

function grow_to!(g :: GeoBucket{K,M}, l ::Int) where {K,M}
    append!(g.tab1,[OrePoly(Vector{K}(undef,2*_d^j),Vector{M}(undef,2*_d^j)) for j in length(g.tab1)+1:l])
    append!(g.tab2,[OrePoly(Vector{K}(undef,2*_d^j),Vector{M}(undef,2*_d^j)) for j in length(g.tab2)+1:l])
    append!(g.active_tab_is_1,[true for j in length(g.active_tab_is_1)+1:l])
    append!(g.indices,[0 for j in length(g.indices)+1:l])
    nothing
end

function rem_lt!(g :: GeoBucket, m :: M) where M 
    indices = g.indices
    for i in length(g):-1:1 
        if indices[i] > 0 && lm(g,i) == m
            indices[i] -= 1 
            return nothing
        end
    end 
    # should never happen 
    return nothing 
end 

function prettyprint(g :: GeoBucket,A :: OreAlg)
    for i in 1:length(g)
        tab = active_tab(g,i)

        if iszero(g,i) 
            println("(0)")
            continue
        end
        print("(",coeff(tab,g.indices[i]),")")
        printmon(mon(tab,g.indices[i]),A)
        for j in g.indices[i]-1:-1:1
            print(" + (", coeff(tab,j),")")
            printmon(mon(tab,j),A)
        end
        println()
    end
end


function normalform(g :: GeoBucket{K,M},A :: OreAlg) where {K,M}
    j = 0 
    i= 1 
    while i <= length(g) # length(g) can increase during the loop
        if !iszero(g,i)
            if j != 0 
                add_to_slice_j!(g,j,i,A)
            end
            if !iszero(g,i)
                j=i
            else 
                j = 0 
            end 
        end
        i += 1
    end
    if j == 0
        return zero(A)
    end

    n  = g.indices[j]

    res = undefOrePoly(n,A)
    tab = active_tab(g,j)
    for i in 1:n
        res[i] = tab[n-i+1]
    end
    g.indices[j] = 0
    return res
end

function lm(g :: GeoBucket, A :: OreAlg)
    while true
        j = 0 
        for i in 1:length(g)
            if !iszero(g,i)
                if j == 0 || lt(order(A), lm(g,j),lm(g,i))
                    j=i 
                elseif lm(g,i) == lm(g,j)
                    tabj = active_tab(g,j)
                    ind = g.indices[j]
                    tabj.coeffs[ind] = add(lc(g,i), coeff(tabj,ind), ctx(A))
                    g.indices[i] -= 1
                    if iszero(tabj.coeffs[g.indices[j]])
                        g.indices[j] -= 1
                        j = -1
                        break
                    end
                end
            end
        end

        if j == 0 
            error("taking lm of a zero polynomial. If unsure you should use lt and check if the coefficient is non-zero")
            return makemon(-1,A)
        elseif j > 0 
        return lm(g,j)
        end
    end
end



function lt(g :: GeoBucket, A :: OreAlg)
    while true
        j = 0 
        for i in 1:length(g)
            if !iszero(g,i)
                if j == 0 || lt(order(A), lm(g,j),lm(g,i))
                    j=i 
                elseif lm(g,i) == lm(g,j)
                    tabj = active_tab(g,j)
                    ind = g.indices[j]
                    tabj.coeffs[ind] = add(lc(g,i), coeff(tabj,ind), ctx(A))
                    g.indices[i] -= 1
                    if iszero(tabj.coeffs[g.indices[j]])
                        g.indices[j] -= 1
                        j = -1
                        break
                    end
                end
            end
        end

        if j == 0 
            return (zero(ctx(A)),makemon(-1,A))
        elseif j > 0 
            return (lc(g,j),lm(g,j))
        end
    end
end


function add_to_slice_j!(g :: GeoBucket,i :: Int, j :: Int, A :: OreAlg)
    if j > length(g)
        grow_to!(g,j)
    end
    res = inactive_tab(g,j)
    tabi = active_tab(g,i)
    tabj = active_tab(g,j)

    # merge 
    ilen = g.indices[i] 
    jlen = g.indices[j] 
    ii = 1
    ij = 1 
    w = 1 # writing index 

    while (ii <= ilen) && (ij <= jlen)
        mi = mon(tabi,ii)
        mj = mon(tabj,ij)
        if lt(order(A),mi,mj)
            res[w] = tabi[ii]
            ii +=1 
            w += 1 
        elseif mi == mj  
            c = add(coeff(tabi,ii), coeff(tabj,ij),ctx(A))
            if !iszero(c,ctx(A))
                res[w] = (c, mi)
                w +=1 
            end
            ii += 1 
            ij += 1 

        else 
            res[w] = tabj[ij]
            ij += 1 
            w += 1 
        end
    end
    while ii <= ilen 
        res[w] = tabi[ii]
        ii += 1 
        w += 1 
    end
    while ij <= jlen 
        res[w] = tabj[ij]
        ij +=1 
        w += 1 
    end
    g.indices[i] = 0
    g.indices[j] = w - 1

    change_active_tab!(g,j)
    if w > _d^j + 1 
        add_to_slice_j!(g,i+1,i+2,A)
    end

    return g
end

function add!(g :: GeoBucket,co ::T, p :: OrePoly,A :: OreAlg) where T
    j =  int_log_upper(_d,length(p))
    if j > length(g)
        grow_to!(g,j)
    end

    res = inactive_tab(g,j)
    tabj = active_tab(g,j)

    # merge 
    jlen = g.indices[j] 
    ip = length(p)
    ij = 1 
    w = 1 # writing index 
    
    while (ip >=1) && (ij <= jlen)
        mp = mon(p,ip)
        mj = mon(tabj,ij)
        if lt(order(A),mp,mj)
            c = mul(co,coeff(p,ip),ctx(A))
            res[w] = (c,mp)
            ip -=1 
            w += 1 
        elseif mp == mj  
            c = add(mul(co,coeff(p,ip),ctx(A)), coeff(tabj,ij),ctx(A))
            if !iszero(c,ctx(A))
                res[w] = (c, mp)
                w +=1 
            end
            ip -= 1 
            ij += 1 

        else 
            res[w] = tabj[ij]
            ij += 1 
            w += 1 
        end
    end
    while ip >= 1 
        c = mul(co,coeff(p,ip),ctx(A))
        mp = mon(p,ip)
        res[w] = (c,mp)
        ip -= 1 
        w += 1 
    end
    while ij <= jlen 
        res[w] = tabj[ij]
        ij +=1 
        w += 1 
    end
    g.indices[j] = w - 1
    change_active_tab!(g,j)
    if w > _d^j + 1 
        add_to_slice_j!(g,j,j+1,A)
    end
    return g
end
    






# add the product c*m*f to g. 
# The function used depends on the choice of representation of the monomials 
function addmul_geobucket!(g :: GeoBucket, c:: K, m :: M, f :: OrePoly, A:: OreAlg; mod_der :: Val{B} = Val(false) ) where {K,M,B}
    if (A.nlv > 0) && locvar_derivative_interaction(m, f, A)
        # Fallback to the generic multiplication when localisation variables
        # interact with derivatives so that we keep the correct chain rule terms.
        current = normalform(g, A)
        prod = mul((c, m), f, A)
        res = add(current, prod, A)
        fill!(g.indices, 0)
        init_GeoBucket!(g, res)
        return g
    end
    if A.varord == "dright" 
        return addmul_geobucket_dright!(g :: GeoBucket, c:: K, m :: M, f :: OrePoly, A:: OreAlg)
    else # A.varord = "dleft" 
        return addmul_geobucket_dleft!(g :: GeoBucket, c:: K, m :: M, f :: OrePoly, A:: OreAlg, mod_der=Val(B))
    end
end

function locvar_derivative_interaction(m :: M, f :: OrePoly, A :: OreAlg) where M
    der_end = A.nrdv + A.npdv
    loc_start = A.nrdv + 2*A.npdv + A.npv + 1
    loc_end = loc_start + A.nlv - 1

    der_end == 0 && loc_start > loc_end && return false

    has_der_m = any(m[i] > 0 for i in 1:der_end)
    has_loc_m = loc_start <= loc_end && any(m[i] > 0 for i in loc_start:loc_end)

    has_der_f = false
    has_loc_f = false
    for i in 1:length(f)
        mo = mon(f, i)
        if !has_der_f
            for j in 1:der_end
                if mo[j] > 0
                    has_der_f = true
                    break
                end
            end
        end
        if !has_loc_f && loc_start <= loc_end
            for j in loc_start:loc_end
                if mo[j] > 0
                    has_loc_f = true
                    break
                end
            end
        end
        has_der_f && has_loc_f && break
    end

    return (has_loc_m && has_der_f) || (has_der_m && has_loc_f)
end

# this new type is just here to create a custom iterator needed in addmul_geobucket 
struct itr_mon_dright{N,E}
    m  :: SVector{N,E}
    nrdv :: Int
    npdv :: Int 
end

function Base.iterate(m::itr_mon_dright{N,E}) where {N,E}
    return (OreMonVE{N,E}(m.m),m.m)
end

function Base.iterate(m::itr_mon_dright{N,E},v :: SVector{N,E}) where {N,E}
    j = m.nrdv+1
    bound = m.nrdv +m.npdv 
    while j <= bound
        if v[j] > 0
            d = gen_mon_dright(j, m,v)
            return (OreMonVE{N,E}(d),d)
        end
        j += 1
    end
    return nothing 
end

# I get type unstability if the svector creation is not embed in a function ??? 
function gen_mon_dright(j :: Int,m::itr_mon_dright{N,E},v :: SVector{N,E}) where {N,E}
    return SVector{N,E}(i < j ? m.m[i] : (i == j) ? v[i] - E(1) : v[i] for i in 1:N)
end

function guess_size_dright(mm :: M, m :: M, g :: OrePoly, A :: OreAlg) where M 
    res = 0 
    for mmm in mons(g)
        for i in A.nrdv+1:A.nrdv+A.npdv 
            if m[i]-mm[i] > mmm[i+A.npdv]
                @goto next
            end
        end
        res += 1
        @label next
    end
    return res
end

function addmul_geobucket_next_term_dright(i_f :: Int,c :: K, m::M, mm :: M, f :: OrePoly, A::OreAlg) where {K,M}

    while i_f >= 1 
        mf = mon(f,i_f)
        for i in A.nrdv+1:A.nrdv+A.npdv 
            if m[i]-mm[i] > mf[i+A.npdv]
                i_f -=1
                @goto next2
            end
        end
        cc = mul(c,coeff(f,i_f),ctx(A))
        for i in A.nrdv+1:A.nrdv+A.npdv
            expmd = Int(m[i])
            expmfx = Int(mf[i + A.npdv]) 
            j = expmd - mm[i]
            tmp = binomial(expmd,j)
            if j != 0 
                tmp = tmp * prod(expmfx-i for i in 0:j-1)
            end
            tmp = convertn(tmp,ctx(A))
            cc = mul(cc,tmp,ctx(A))
        end
        N = nvars(A)
        mmm = SVector{N,Int16}(i < A.nrdv + A.npdv + 1 || i > A.nrdv + 2*A.npdv ? mf[i] : mf[i] - (m[i-A.npdv] - mm[i-A.npdv]) for i in 1:N)
        mmm = OreMonVE{N,Int16}(mmm)

        return (cc,mm*mmm,i_f-1) 
        @label next2
    end
    return (zero(ctx(A)),makemon(-1,A), -1)
end

function addmul_geobucket_dright!(g :: GeoBucket, c:: K, m :: M, f :: OrePoly, A:: OreAlg ) where {K,M}
    for (l,mm) in enumerate(itr_mon_dright(m.exp,A.nrdv,A.npdv)) 
        siz = guess_size_dright(mm,m,f,A)
        i =  int_log_upper(_d,siz)
        if i > length(g)
            grow_to!(g,i)
        end

        res = inactive_tab(g,i)
        tab = active_tab(g,i)

        # starting to merge in res the array tab with another array computed on the fly 
        (co,mo,i_f) = addmul_geobucket_next_term_dright(length(f),c,m, mm, f, A) # terms computed on the fly
        if i_f < 0 
            continue 
        end

        leng = g.indices[i] 
        ig = 1
        mg = mon(tab,ig)
        w = 1

        while ig <= leng && i_f >= 0
            mg = mon(tab,ig)
            if lt(order(A),mg,mo)
                res[w] = tab[ig]
                ig += 1 
                w +=1
            elseif mo == mg 
                cc = add(co, coeff(tab,ig),ctx(A))
                if !iszero(cc,ctx(A))
                    res[w] = (cc,mo)
                    w += 1 
                end
                (co, mo, i_f) = addmul_geobucket_next_term_dright(i_f,c,m, mm, f, A)
                ig += 1
            else 
                res[w] = (co,mo)
                (co, mo, i_f) = addmul_geobucket_next_term_dright(i_f,c,m, mm, f, A)
                w += 1
            end
        end
        while ig <= leng
            res[w] = tab[ig]
            ig += 1 
            w += 1 
        end

        while i_f >= 0
            res[w] = (co,mo)
            (co, mo, i_f) = addmul_geobucket_next_term_dright(i_f,c,m, mm, f, A)
            w += 1
        end


        g.indices[i] = w-1
        change_active_tab!(g,i)
        if w > _d^i + 1
            add_to_slice_j!(g,i,i+1,A)
        end
    end
    return g
end


# same with dleft


# this new type is just here to create a custom iterator needed in addmul_geobucket 
struct itr_mon_dleft{N,E}
    m  :: SVector{N,E}
    nrdv :: Int
    npdv :: Int 
end

function Base.iterate(m::itr_mon_dleft{N,E}) where {N,E}
    return (OreMonVE{N,E}(m.m),m.m)
end

function Base.iterate(m::itr_mon_dleft{N,E},v :: SVector{N,E}) where {N,E}
    j = m.nrdv+m.npdv+1
    bound = m.nrdv + 2*m.npdv 
    while j <= bound
        if v[j] > 0
            d = gen_mon_dleft(j, m,v)
            return (OreMonVE{N,E}(d),d)
        end
        j += 1
    end
    return nothing 
end

# I get type unstability if the svector creation is not embed in a function ??? 
function gen_mon_dleft(j :: Int,m::itr_mon_dleft{N,E},v :: SVector{N,E}) where {N,E}
    return SVector{N,E}(i < j ? m.m[i] : (i == j) ? v[i] - E(1) : v[i] for i in 1:N)
end

function guess_size_dleft(mm :: M, m :: M, g :: OrePoly, A :: OreAlg) where M 
    res = 0 
    for mmm in mons(g)
        for i in A.nrdv+A.npdv+1:A.nrdv+2*A.npdv 
            if m[i]-mm[i] > mmm[i-A.npdv]
                @goto next
            end
        end
        res += 1
        @label next
    end
    return res
end

function addmul_geobucket_next_term_dleft(i_f :: Int,c :: K, m::M, mm :: M, f :: OrePoly, A::OreAlg) where {K,M}

    while i_f >= 1 
        mf = mon(f,i_f)
        for i in A.nrdv+A.npdv+1:A.nrdv+2*A.npdv 
            if m[i]-mm[i] > mf[i-A.npdv]
                i_f -=1
                @goto next2
            end
        end
        cc = mul(c,coeff(f,i_f),ctx(A))
        for i in A.nrdv+A.npdv+1:A.nrdv+2*A.npdv
            expmx = Int(m[i])
            expmfd = Int(mf[i - A.npdv]) 
            j = expmx - mm[i]
            tmp = binomial(expmx,j)
            if j != 0 
                tmp = tmp * prod(expmfd-i for i in 0:j-1)
            end
            tmp = convertn(tmp,ctx(A))
            if isodd(j)
                tmp = opp(tmp,ctx(A))
            end
            cc = mul(cc,tmp,ctx(A))
        end
        N = nvars(A)
        mmm = SVector{N,Int16}(i < A.nrdv + 1 || i > A.nrdv + A.npdv ? mf[i] : mf[i] - (m[i+A.npdv] - mm[i+A.npdv]) for i in 1:N)
        mmm = OreMonVE{N,Int16}(mmm)

        return (cc,mm*mmm,i_f-1) 
        @label next2
    end
    return (zero(ctx(A)),makemon(-1,A), -1)
end


function addmul_geobucket_dleft!(g :: GeoBucket, c:: K, m :: M, f :: OrePoly, A:: OreAlg; mod_der :: Val{B} = Val(false)) where {K,M,B}
    for (l,mm) in enumerate(itr_mon_dleft(m.exp,A.nrdv,A.npdv)) 
        siz = guess_size_dleft(mm,m,f,A)
        i =  int_log_upper(_d,siz)
        if i > length(g)
            grow_to!(g,i)
        end

        res = inactive_tab(g,i)
        tab = active_tab(g,i)

        # starting to merge in res the array tab with another array computed on the fly 
        (co,mo,i_f) = addmul_geobucket_next_term_dleft(length(f),c,m, mm, f, A) # terms computed on the fly
        if i_f < 0 
            continue 
        end

        leng = g.indices[i] 
        ig = 1
        mg = mon(tab,ig)
        w = 1

        while ig <= leng && i_f >= 0
            mg = mon(tab,ig)
            if lt(order(A),mg,mo)
                res[w] = tab[ig]
                ig += 1 
                w +=1
            elseif mo == mg 
                cc = add(co, coeff(tab,ig),ctx(A))
                if !iszero(cc,ctx(A))
                    res[w] = (cc,mo)
                    w += 1 
                end
                (co, mo, i_f) = addmul_geobucket_next_term_dleft(i_f,c,m, mm, f, A)
                ig += 1
            else 
                res[w] = (co,mo)
                (co, mo, i_f) = addmul_geobucket_next_term_dleft(i_f,c,m, mm, f, A)
                w += 1
            end
        end
        while ig <= leng
            res[w] = tab[ig]
            ig += 1 
            w += 1 
        end

        while i_f >= 0
            res[w] = (co,mo)
            (co, mo, i_f) = addmul_geobucket_next_term_dleft(i_f,c,m, mm, f, A)
            w += 1
        end

        g.indices[i] = w-1
        change_active_tab!(g,i)
        B && mod_derivatives!(g,i,A)
        w = g.indices[i]
        if w > _d^i
            add_to_slice_j!(g,i,i+1,A)
        end
    end
    return g
end


###

    
function int_log_upper(d::Int,a::Int)
    res = d 
    pow = 1
    while res <= a 
        res = res * d 
        pow += 1 
    end 
    return pow 
end

function nozero(g :: GeoBucket, A :: OreAlg)
    for i in 1:length(g)
        tab = active_tab(g,i)
        for j in 1:g.indices[i]
            iszero(coeff(tab,j),ctx(A)) && return false
        end
    end
    return true
end


# It assumes that every pivots have been found
function reducebasis!(geob :: GeoBucket, tmp_poly :: ReuseOrePoly, f :: Vector{OrePoly{K,M}}, A :: Alg) where {K,M, Alg <:OreAlg}
    i = 1

    lenf = length(f)
    
    while i <= lenf
        makemonic!(f[i],A)
        tmp = f[lenf]
        f[lenf] = f[i]
        f[i] = tmp
        gb = @view f[1:lenf-1]
        f[lenf] = div!(geob, tmp_poly,f[lenf], gb, A)
        if length(f[lenf]) == 0 
            lenf = lenf -1 
            pop!(f)
        else
            tmp = f[lenf]
            f[lenf] = f[i]
            f[i] = tmp
            i = i + 1
        end
    end

    sort!(f, lt = (x,y) -> lt(order(A),x[1][2], y[1][2]), rev = true)
end

function div!(geob :: GeoBucket ,tmp_poly :: ReuseOrePoly,f :: OrePoly, gb :: Vector{OrePoly{K,M}}, A :: OreAlg) where {K,M}
    geob = init_GeoBucket!(geob,f)
    while true 
        div = false
        lco,lmon = lt(geob,A)
        if iszero(lco)
            break
        end
        for i in length(gb):-1:1
            lmon2 = lm(gb[i])
            if iscompatible(lmon2, lmon,A) && divide(lmon2, lmon,A)
                div = true
                geob = reduce_geob!(geob, (lco,lmon), gb[i], A,Val(false))
                break
            end
        end
        if !div
            push!(tmp_poly,lco,lmon)
            rem_lt!(geob,lmon)
        end
    end
    # res = normalform(geob,A)
    return copy_to_OrePoly!(tmp_poly,A)
end


function reduce_geob!(g :: GeoBucket, t:: Tuple{K,M}, f :: OrePoly{K,M}, A :: OreAlg, :: Val{B})  where {K, M, B}
    (c,m) = t # lc and lm of g   
    themon = m/lm(f) 
    thectx = ctx(A)
    thecoeff = opp(mul(c,inv(lc(f),thectx),thectx),thectx)

    addmul_geobucket!(g, thecoeff, themon, f, A, mod_der = Val(B))
    return g 
end
