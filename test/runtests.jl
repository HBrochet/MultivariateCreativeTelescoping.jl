using MultivariateCreativeTelescoping
using Test

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

    #testing commute 
    A = OreAlg(order = "lex d1 d2 dx x",ratdiffvars=(["t1","t2"],["d1","d2"]),poldiffvars=(["x"],["dx"]),char = 373587883)
    p = parse_OrePoly("t1*x",A)
    q = parse_OrePoly("d2",A)
    @test commute(p,q,A)

    p = parse_OrePoly("t1*x",A)
    q = parse_OrePoly("d1",A)
    @test !commute(p,q,A)
end