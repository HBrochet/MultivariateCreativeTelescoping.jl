struct GeoBucket{K, M <: AbsOreMonomial} <: AbsOrePolynomial{K, M}
    tab1 :: Vector{OrePoly{K,M}}
    tab2 :: Vector{OrePoly{K,M}}
    active_tab_is_1 :: Vector{Bool} # contains the value true (resp. false) is tab1[i] is active (resp. tab2[i])
    indices :: Vector{Int} # index of the lm in the ith active tab 
end

# following Yan's paper we create an array tab s.t length(tab[i]) <= 4*4^i 
# that store a polynomial of size at most 4^i
# It is stored in the first 4^i indices if indices[i] <= 4^i and between the indices 2*4^i +1:end 
# otherwise. The remaining memory is preallocated for addition of polynomials 
# !!! Monomials are stored in increasing order in GeoBuckets

function GeoBucket(p :: OrePoly{K,M}) where {K, M<: AbsOreMonomial}
    lenp = length(p)
    l = int_log_upper(4,lenp)
    tab1 = [OrePoly(Vector{K}(undef,2*4^i),Vector{M}(undef,2*4^i)) for i in 1:l]
    tab2 = [OrePoly(Vector{K}(undef,2*4^i),Vector{M}(undef,2*4^i)) for i in 1:l]
    active_tab_is_1 = [true for j in 1:l]
    indices = [0 for j in 1:l]
    for j in 1:lenp
        tab1[l][lenp-j+1] = p[j]
    end
    indices[l] = lenp
    return GeoBucket{K,M}(tab1,tab2,active_tab_is_1,indices) 
end

tab1_is_active(g :: GeoBucket, i :: Int) = g.active_tab_is_1[i]
Base.length(g :: GeoBucket) = length(g.tab1)
Base.iszero(g :: GeoBucket, i ::Int) = g.indices[i] == 0
 
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
    append!(g.tab1,[OrePoly(Vector{K}(undef,2*4^j),Vector{M}(undef,2*4^j)) for j in length(g.tab1)+1:l])
    append!(g.tab2,[OrePoly(Vector{K}(undef,2*4^j),Vector{M}(undef,2*4^j)) for j in length(g.tab2)+1:l])
    append!(g.active_tab_is_1,[true for j in length(g.active_tab_is_1)+1:l])
    append!(g.indices,[0 for j in length(g.indices)+1:l])
    nothing
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


# add the product c*m*f to g  
function addmul_geobucket!(g :: GeoBucket, c:: K, m :: M, f :: OrePoly, A:: OreAlg ) where {K,M}
    supp =  supp_mons(m,A) 
    # supp is a set of monomials 
    # mf is the sum of multiples OrePoly (f_n)_{n\in supp}
    # where f_n has a support in n*supp(f)  (this multiplication is commutative, ie exponents are added)
    # and size sizes[i] if n = supp[i]
    # these polynomials f_n are added on the fly to g 

    sizes = guess_sizes(supp,m,f,A) 

    ### lets go 
    for (l,mm) in enumerate(supp) 
        i =  int_log_upper(4,sizes[l])
        if i > length(g)
            grow_to!(g,i)
        end

        res = inactive_tab(g,i)
        tab = active_tab(g,i)

        # starting to merge in res the array tab with another array computed on the fly 
        (co,mo,i_f) = addmul_geobucket_next_term(length(f),c,m, mm, f, A) # terms computed on the fly
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
                (co, mo, i_f) = addmul_geobucket_next_term(i_f,c,m, mm, f, A)
                ig += 1
            else 
                res[w] = (co,mo)
                (co, mo, i_f) = addmul_geobucket_next_term(i_f,c,m, mm, f, A)
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
            (co, mo, i_f) = addmul_geobucket_next_term(i_f,c,m, mm, f, A)
            w += 1
        end


        g.indices[i] = w-1
        change_active_tab!(g,i)
        if w > 4^i + 1
            add_to_slice_j!(g,i,i+1,A)
        end
    end
    return g
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
    if w > 4^j + 1 
        add_to_slice_j!(g,i+1,i+2,A)
    end

    return g
end

function addmul_geobucket_next_term(i_f :: Int,c :: K, m::M, mm :: M, f :: OrePoly, A::OreAlg) where {K,M}

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
            expmd = K(m[i])
            expmfx = K(mf[i + A.npdv]) 
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
        # @assert !iszero(cc)

        return (cc,mm*mmm,i_f-1) 
        @label next2
    end
    return (zero(ctx(A)),makemon(-1,A), -1)
end



function guess_sizes(supp :: Vector{M}, m :: M, g :: OrePoly,A :: OreAlg) where M 
    res = [0 for i in 1:length(supp)]
    for (j,mm) in enumerate(supp) 
        for mmm in mons(g)
            for i in A.nrdv+1:A.nrdv+A.npdv 
                if m[i]-mm[i] > mmm[i+A.npdv]
                    @goto next 
                end
            end
            res[j] += 1
            @label next
        end
    end
    return res
end

function supp_mons(m :: M, A :: OreAlg) where M
    return supp_mons_rec(m, A.nrdv, A)
end

function supp_mons_rec(m :: M, l :: Int, A :: OreAlg) where M
    res = M[] 
    if l+1 == A.nrdv+A.npdv+1 
        return [m] 
    end
    mm = makemon(l+1,A)
    for i in 0:m[l+1] 
        append!(res,supp_mons_rec(m/mm^i,l+1,A))
    end
    return res  
end


    
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