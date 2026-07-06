using StaticArrays
using Nemo: polynomial_ring, degree, Native, fpPolyRingElem

@testset "Cauchy interpolation" begin 
    A = OreAlg(order = "lex dt x dx",ratdiffvars=(["t"],["dt"]), poldiffvars=(["x"],["dx"]),char=primes[1])

    p = parse_OrePoly("t^10//(t+1) + t^2*x",A)

    function foo(A::OreAlg,p :: OrePoly)
        cs = [mul(c,c,ctx(A)) for c in coeffs(p)]
        return OrePoly(cs,deepcopy(mons(p)))
    end
    
    res_default = compute_with_cauchy_interpolation(foo,A,p)
    res_same_den = compute_with_cauchy_interpolation(foo,A,p; param = ci_param(same_den = Val(true)))
    @test res_default == parse_OrePoly("(t^4)*x*1 + (t^2//(t^2 + 2*t + 1))*1",A)
    @test res_same_den == res_default
end

@testset "Cauchy interpolation (kreg-se-ll-4)" begin
    Rpt = OreAlg(order = "lex dt > grevlex p1 p2 p3 p4 > grevlex d1 d2 d3 d4",
                 ratdiffvars = (["t"], ["dt"]),
                 poldiffvars = (["p1", "p2", "p3", "p4"], ["d1", "d2", "d3", "d4"])
                )
    Gstring = "[0+(p1-1/6*t*p1^3+t*p1+1/2*t*p1*p2-1/3*t*p3)*1+(-1/2*t*p1)*d1^2+(-t)*d1*d2+(1/2*t*p1^2-t-1/2*t*p2)*d1+(t*p1)*d2+(t)*d3+(1/6*t)*d1^3, 0+(p2-1/2*t*p1^2+t+1/2*t*p2)*1+(-1/2*t)*d1^2+(t*p1)*d1+(t)*d2, 0+(t)*d1+(-t*p1+p3)*1, 0+(p4-t)*1, 0+(-1-1/3*p1*p3+1/4*p1^2*p2-1/24*p1^4+1/2*p1^2-1/2*p2+1/4*p4-1/8*p2^2)*1+(1)*dt+(1/2-1/4*p1^2+1/4*p2)*d1^2+(-p1)*d1*d2+(-p1-1/2*p1*p2+1/6*p1^3+1/3*p3)*d1+(-1+1/2*p1^2-1/2*p2)*d2+(p1)*d3+(1)*d4+(1/6*p1)*d1^3+(1/2)*d1^2*d2+(-1/2)*d2^2+(-1)*d1*d3+(-1/24)*d1^4]"
    G = parse_vector_OrePoly(Gstring, Rpt)
    B = f4(G, Rpt)
    res = MCT(one(Rpt), B, Rpt)
    res_expected = parse_OrePoly("(16*t^11 + 64*t^10 + 16*t^9 - 128*t^8 + 128*t^7 + 224*t^6 - 640*t^5 - 192*t^4 + 768*t^3 - 256*t^2)*dt^2*1 + (-4*t^13 - 16*t^12 + 64*t^10 + 40*t^9 + 144*t^8 + 880*t^7 + 1392*t^6 + 192*t^5 - 800*t^4 + 1344*t^3 + 960*t^2 - 1664*t + 384)*dt*1 + (-t^14 - 4*t^13 - 4*t^12 - 4*t^11 - 24*t^10 - 24*t^9 + 12*t^8 - 32*t^7 - 48*t^6 + 64*t^5 - 16*t^4)*1", Rpt)
    @test res == res_expected
end

@testset "half_gcd (flint vs julia)" begin
    p = primes[1]
    R, x = polynomial_ring(Native.GF(p), "x")

    a = x^8 + 3*x^5 + 2*x + 1
    b = x^6 + 2*x^3 + x + 1

    Mj = half_gcd_julia(a, b)
    Mf = half_gcd(a, b)

    v = SVector{2,fpPolyRingElem}(a, b)
    rj = Mj * v
    rf = Mf * v

    @test rj == rf
    @test Nemo.degree(rf[2]) < div(Nemo.degree(a), 2) + isodd(Nemo.degree(a))
end
