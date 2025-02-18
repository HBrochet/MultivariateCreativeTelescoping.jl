using MultivariateCreativeTelescoping
using Test

@testset "unused type" begin
    @test isempty(detect_unbound_args(MultivariateCreativeTelescoping)) 
end

@testset "Algebra definition, parsing, add, mul" begin
    A = OreAlg(char = 373587883,order = "lex x1 x2 d1 d2",poldiffvars=(["x1","x2"],["d1","d2"]))    
    pol1 = parse_OrePoly("x1*d1",A)
    pol2 = parse_OrePoly("d1*x1 - 1",A)
    @test pol1 == pol2

    pol3 = parse_OrePoly("x1^2*x2^2*d1^2*d2^2",A)

    s="1*d1^2*d2^2*x1^2*x2^2 + 373587879*d1^2*d2*x1^2*x2 + 2*d1^2*x1^2 + 373587879d1*d2^2*x1*x2^2 + 16*d1*d2*x1*x2 + 373587875*d1*x1 + 2*d2^2*x2^2 + 373587875d2*x2 + 4"
    pol4 = parse_OrePoly(s,A)

    @test pol3 == pol4
    # testing ratdiffvars 
    A = OreAlg(order="lex d1 d2 dx x",ratdiffvars=(["t1","t2"],["d1","d2"]),poldiffvars=(["x"],["dx"]),char = 373587883)
    pol4b = parse_OrePoly("2*x*d1/t1*dx",A)
    pol4c = parse_OrePoly("2/t1*d1*dx*x - 2/t1*d1 - 2/t1^2*dx*x + 2/t1^2",A)

    @test pol4b == pol4c
    #testing loc vars
    A = OreAlg(order = "lex Tx Ty x y dx dy",poldiffvars=(["x","y"],["dx","dy"]),locvars=(["Tx","Ty"],["x^2","y^2"]),char = 373587883)   
    pol5 = parse_OrePoly("Tx^2*dx",A)
    pol6 = parse_OrePoly("4*x*Tx^3 + 1*dx*Tx^2",A)
    @test pol5 == pol6 

    pol7 = parse_OrePoly("Tx*dx^2",A)
    pol8 = parse_OrePoly("8*x^2*Tx^3 + 4*dx*x*Tx^2 - 2*Tx^2 + 1*dx^2*Tx",A)
    @test pol7 == pol8

    #testing ratvars 
    A = OreAlg(order = "lex x y dx dy",ratvars=["s1","s2"],poldiffvars=(["x","y"],["dx","dy"]))   
    pol9 = parse_OrePoly("x*s1*(1-s2)",A)
    @test pol9 isa OrePoly

    #testing block order 
    A = OreAlg(order = "lex dz z > lex x > grevlex y > grevlex dx dy",poldiffvars=(["x","y","z"],["dx","dy","dz"]))    
    q = parse_OrePoly("x*x + y*y + z*z + dz*z+dx*dy + dy*dy + dx*y ",A) 


    #testing commute 
    A = OreAlg(order = "lex d1 d2 dx x",ratdiffvars=(["t1","t2"],["d1","d2"]),poldiffvars=(["x"],["dx"]),char = 373587883)
    p = parse_OrePoly("t1*x",A)
    q = parse_OrePoly("d2",A)
    @test commute(p,q,A)

    p = parse_OrePoly("t1*x",A)
    q = parse_OrePoly("d1",A)
    @test !commute(p,q,A)
end


@testset "F4, F5 and saturation" begin 
    A = OreAlg(order = "lex T > grevlex x y > grevlex dx dy",poldiffvars=(["x","y"],["dx","dy"]),locvars=(["T"],["x^2-y^3"]),char = 373587883,nomul=["T"])
    p = parse_OrePoly("x^2-y^3",A)
    g = ann_inv_pol(p,A)
    @test check_ann(g,p,A)

    #same computation with different coefficient type 
    A = OreAlg(order = "lex T > grevlex x y > grevlex dx dy dt",ratdiffvars=(["t"],["dt"]),poldiffvars=(["x","y"],["dx","dy"]),locvars=(["T"],["x^2-y^3"]),char = 373587883,nomul=["T"])
    p = parse_OrePoly("x^2-y^3",A)
    g = ann_inv_pol(p,A)
    #@test check_ann(g,p,A) # todo: write a check_ann function when coeff = "ratfunmodp"


    A = OreAlg(order = "lex T > grevlex x y z > grevlex dx dy dz",poldiffvars=(["x","y","z"],["dx","dy","dz"]),locvars=(["T"],["x*y^2*z^2+x*y^2+x*z^2-y*z+x"]),char = 373587883,nomul=["T"])
    p = parse_OrePoly("x*y^2*z^2+x*y^2+x*z^2-y*z+x",A)
    g = ann_inv_pol(p,A)
    @test check_ann(g,p,A)
    #"x*y^2*z^2+x*y^2+x*z^2-y*z+x"


    A = OreAlg(order="lex dx da",ratdiffvars=(["a","x"],["da","dx"]),char = 373587883)
    p1 = parse_OrePoly("a*da - x*dx + 1",A)
    p2 = parse_OrePoly("x*dx^2-dx+a^2*x",A)
    g = [p1,p2]
    gb1 = f4(g,A)
    sgb = sigbasis(g,A)
    gb2 = sgbtogb(sgb,A)
    @test gb1 == gb2
end

@testset "Weyl closure and MCT" begin
    A = OreAlg(order = "grevlex y > grevlex dx dy",ratdiffvars=(["x"],["dx"]),poldiffvars=(["y"],["dy"]),char = primes[1],varord = "dright")
    p = parse_OrePoly("x^2-y^3",A)
    gens = [parse_OrePoly("dx*(x^2-y^3)",A),parse_OrePoly("dy*(x^2-y^3)",A)]
    init = weyl_closure_init(A)
    g = weyl_closure(gens,A,init)
    @test isholonomic(g,A)

    A = OreAlg(order = "grevlex x y > grevlex dx dy",poldiffvars=(["x","y"],["dx","dy"]))
    p = parse_OrePoly("x^2-y^3",A)
    gens = [parse_OrePoly("dx*(x^2-y^3)",A),parse_OrePoly("dy*(x^2-y^3)",A)]
    init = weyl_closure_init(A)
    g = weyl_closure(gens,A,init)
    @test isholonomic(g,A)

    # t
    A = OreAlg(order = "lex dt x dx",ratdiffvars=(["t"],["dt"]), poldiffvars=(["x"],["dx"]))
    ann = [parse_OrePoly("dt*(x-t)",A),parse_OrePoly("dx*(x-t)",A)]
    spol = parse_OrePoly("x",A)
    init = weyl_closure_init(A)
    gb = weyl_closure(ann,A,init)
    res = MCT(spol, gb, A)
    res2 = parse_OrePoly("t*dt-1",A)
    @test res == res2


    # # SSW3
    A = OreAlg(order = "lex dt > grevlex x y > grevlex dx dy",ratdiffvars=(["t"],["dt"]),poldiffvars = (["x","y"],["dx","dy"]))
    s = "(t*x^2*y^2+t*x*y^2+t*x^2+t*y^2+t*x-x*y+t)"
    ann = [parse_OrePoly("dt*"*s,A),parse_OrePoly("dx*"*s,A),parse_OrePoly("dy*"*s,A)]
    spol = parse_OrePoly("(-x^2*y^2+x^2+y^2-1)",A)
    init = weyl_closure_init(A)
    gb = weyl_closure(ann,A,init)
    res = MCT(spol, gb, A)


    
    # res2 = parse_OrePoly("t*dt-1",A)
    # @test res == res2

    
    
end


@testset begin 
    A = OreAlg(order = "lex dt x dx",ratdiffvars=(["t"],["dt"]), poldiffvars=(["x"],["dx"]),char=primes[1])

    p = parse_OrePoly("t",A)

    function foo(A::OreAlg,p :: OrePoly)
        cs = [mul(c,c,ctx(A)) for c in coeffs(p)]
        return OrePoly(cs,deepcopy(mons(p)))
    end
    
    res = compute_with_cauchy_interpolation(foo,A,p)
    @test res == parse_OrePoly("t^2",A)
end

### test for mri 


# global ctr = 0

# A = OreAlg(order = "grevlex x y dx dy",poldiffvars=(["x","y"],["dx","dy"]),ratvars=["s","t"],char=primes[1])

# p = parse_OrePoly("s^2*(t+1)/(s+1)/(t+2)*y^2 + (t^2+1)*(s+1)/(s^4*t +1)",A)

# q = parse_OrePoly("s^2*(t+1)/(s+1)/(t+2)*y^2 + (t^2+1)*(s+1)/(s^4*t)",A)

# q0 = parse_OrePoly("s^2*(t+1)",A)
# q0b = parse_OrePoly("t^2*(s+1)",A)


# q1 = parse_OrePoly("s^2*(t+1)/t^7",A)

# r = parse_OrePoly("1/s",A)

# s = parse_OrePoly("(s^2 + t^2) /(s^3 - t^2*s)",A)

# t = parse_OrePoly("((s+1)^2 + (t+1)^2)/((s+1)^3 - (t+1)^2*(s+1))",A) # 31

# u = parse_OrePoly("t*(s+1)^10/(t+1)^15",A) #378

# v = parse_OrePoly("1/(s^2*t^2)",A) # 23 ?? 

# function test(A,p)
#     global ctr += 1
#     return p
# end

# res = multivariate_rational_interpolation(test,A,q0)



## mri coupled with a GB calculation 

# using MultivariateCreativeTelescoping
# s = "[-e*m1^2*x1^2-e*m1^2*x1*x2+e*m2^2*x1^2+e*m2^2*x1*x2-e*t*x1^2+e*t*x1*x2+m1^2*x1^2+m1^2*x1*x2+m2^2*x1^2+3*m2^2*x1*x2+2*m2^2*x2^2-t*x1^2-t*x1*x2+(m1^2*x1^3+2*m1^2*x1^2*x2+m1^2*x1*x2^2+m2^2*x1^2*x2+2*m2^2*x1*x2^2+m2^2*x2^3-t*x1^2*x2-t*x1*x2^2)*dx2, e*m1^2*x1*x2+e*m1^2*x2^2-e*m2^2*x1*x2-e*m2^2*x2^2+e*t*x1*x2-e*t*x2^2+2*m1^2*x1^2+3*m1^2*x1*x2+m1^2*x2^2+m2^2*x1*x2+m2^2*x2^2-t*x1*x2-t*x2^2+(m1^2*x1^3+2*m1^2*x1^2*x2+m1^2*x1*x2^2+m2^2*x1^2*x2+2*m2^2*x1*x2^2+m2^2*x2^3-t*x1^2*x2-t*x1*x2^2)*dx1, -x1*x2*e-x1*x2+(m1^2*x1^2+m1^2*x1*x2+m2^2*x1*x2+m2^2*x2^2-t*x1*x2)*dt]"
# s2 = "(x2+x1)*(m1^2*x1^2+m1^2*x1*x2+m2^2*x1*x2+m2^2*x2^2-t*x1*x2)"

# A = OreAlg(char=primes[1],order = "lex T dt > grevlex x1 x2 > grevlex dx1 dx2",ratdiffvars =(["t"],["dt"]),poldiffvars=(["x1","x2"],["dx1","dx2"]),locvars=(["T"],[s2]),ratvars=["e","m1","m2"],nomul=["dt","T"])

# gens0 = parse_vector_OrePoly(s,A)

# gens = deepcopy(gens0)
# T = parse_OrePoly("T",A)
# l = length(gens)
# for i in 1:l 
#     push!(gens, mul(T,gens[i],A))
# end
# push!(gens, parse_OrePoly(s2*"T-1",A))
# push!(gens, parse_OrePoly(s2*"T^2-T",A))
# push!(gens, parse_OrePoly("dt*("*s2*"*T-1)",A))

# param = f5_param(stophol = Val(true))

# gb = f5(gens,A,param=param)

# gb2 = multivariate_rational_interpolation(f5_mri,A,gens,param)
# gb2 = unflatten(gb2,A)