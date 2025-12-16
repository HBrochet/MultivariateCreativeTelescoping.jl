using Test
using Nemo: QQ
using MultivariateCreativeTelescoping

@testset "OrePolyAddMul normalize!" begin
    A = OreAlg(order = "lex x dx", poldiffvars = (["x"], ["dx"]))
    dx = makemon(A.strvar_to_indexp["dx"], A)
    x = makemon(A.strvar_to_indexp["x"], A)

    poly = OrePoly([QQ(1), QQ(2), QQ(-2)], [dx, x, x])
    normal_poly = @inferred normalize!(poly, A)
    expected = parse_OrePoly("dx", A)
    @test normal_poly == expected
    @test normal_poly === poly

    zero_poly = OrePoly([QQ(1), -QQ(1)], [x, x])
    @inferred normalize!(zero_poly, A)
    @test iszero(zero_poly)
end

@testset "OrePolyAddMul add variants" begin
    A = OreAlg(order = "lex x dx", poldiffvars = (["x"], ["dx"]))
    p = parse_OrePoly("x + dx", A)
    q = parse_OrePoly("2*x - dx", A)
    expected = parse_OrePoly("3*x", A)

    res_mut = @inferred add!(copy(p), q, A)
    @test res_mut == expected

    res_add2 = @inferred add2!(copy(p), q, A)
    @test res_add2 == expected

    res_add2_no = @inferred add2!(copy(p), q, A; normal = false)
    @test length(res_add2_no) == length(p) + length(q)

    res_nonmut = @inferred add(p, q, A)
    @test res_nonmut == expected
    @test p == parse_OrePoly("x + dx", A)

    res_add2_nonmut = @inferred add2(p, q, A)
    @test res_add2_nonmut == expected

    vec = [parse_OrePoly("x", A), parse_OrePoly("dx", A), parse_OrePoly("-dx", A)]
    res_vec = @inferred add(vec, A)
    @test res_vec == parse_OrePoly("x", A)
end

@testset "Geobucket respects localisation variable" begin
    A = OreAlg(order = "lex T > grevlex x > grevlex dx",
               poldiffvars = (["x"], ["dx"]),
               locvars = (["T"], ["x^2+1"]),
               nomul = ["T"],
               char = 373587883)

    dx_pol = parse_OrePoly("dx", A)
    loc_rel = sub!(mul(parse_OrePoly("x^2+1", A), parse_OrePoly("T", A), A), one(A), A)

    lcm_mon = lcm(lm(dx_pol), lm(loc_rel))
    geob = GeoBucket(zero(A))
    addmul_geobucket!(geob, lc(loc_rel), lcm_mon / lm(dx_pol), dx_pol, A)
    addmul_geobucket!(geob, opp(lc(dx_pol), ctx(A)), lcm_mon / lm(loc_rel), loc_rel, A)

    res_geob = normalform(geob, A)
    res_std = sub(mul(lcm_mon / lm(dx_pol), dx_pol, A),
                  mul(lcm_mon / lm(loc_rel), loc_rel, A),
                  A)

    @test res_geob == res_std
end
