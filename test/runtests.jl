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
    A = OreAlg(order = "grevlex y > grevlex dx dy",ratdiffvars=(["x"],["dx"]),poldiffvars=(["y"],["dy"]))
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

    # # # SSW3
    # A = OreAlg(order = "lex dt > grevlex x y > grevlex dx dy",ratdiffvars=(["t"],["dt"]),poldiffvars = (["x","y"],["dx","dy"]))
    # s = "(t*x^2*y^2+t*x*y^2+t*x^2+t*y^2+t*x-x*y+t)"
    # ann = [parse_OrePoly("dt*"*s,A),parse_OrePoly("dx*"*s,A),parse_OrePoly("dy*"*s,A)]
    # spol = parse_OrePoly("(-x^2*y^2+x^2+y^2-1)",A)
    # init = weyl_closure_init(A)
    # gb = weyl_closure(ann,A,init)
    # res = MCT(spol, gb, A)


    
    # res2 = parse_OrePoly("t*dt-1",A)
    # @test res == res2

    
    
end


