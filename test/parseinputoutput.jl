using Test
using Nemo: QQ, QQFieldElem
using MultivariateCreativeTelescoping

@testset "ParseInputOutput parsing" begin
    A = OreAlg(order = "lex x dx", poldiffvars = (["x"], ["dx"]))

    zero_poly = @inferred parse_OrePoly("0", A)
    @test iszero(zero_poly)

    const_poly = @inferred parse_OrePoly("3", A)
    @test coeff(const_poly, 1) == QQ(3)
    @test mon(const_poly, 1) == makemon(-1, A)

    x_poly = @inferred parse_OrePoly("x", A)
    @test coeff(x_poly, 1) == QQ(1)
    @test mon(x_poly, 1) == makemon(A.strvar_to_indexp["x"], A)

    expr_poly = @inferred parse_OrePoly("x + 2*dx - 3", A)
    expected_expr = OrePoly(
        [QQ(1), QQ(2), QQ(-3)],
        [makemon(A.strvar_to_indexp["x"], A),
         makemon(A.strvar_to_indexp["dx"], A),
         makemon(-1, A)],
    )
    normalize!(expected_expr, A)
    @test expr_poly == expected_expr

    neg_poly = @inferred parse_OrePoly("-x", A)
    @test coeff(neg_poly, 1) == -QQ(1)

    div_poly = @inferred parse_OrePoly("x/2", A)
    @test coeff(div_poly, 1) * QQ(2) == QQ(1)

    big_poly = @inferred parse_OrePoly("big\"12345678901234567890\"", A)
    @test coeff(big_poly, 1) == QQ(parse(BigInt,"12345678901234567890"))
    @test mon(big_poly, 1) == makemon(-1, A)

    vec = @inferred parse_vector_OrePoly("[x, dx, 1]", A)
    @test vec isa Vector{OrePoly{QQFieldElem, eltype_mo(A)}}
    @test length(vec) == 3
    @test vec[1] == x_poly
    @test vec[2] == parse_OrePoly("dx", A)
    @test coeff(vec[3], 1) == QQ(1)
end

@testset "ParseInputOutput string roundtrip" begin
    A = OreAlg(order = "lex x dx", poldiffvars = (["x"], ["dx"]))
    poly = parse_OrePoly("x + dx - 2", A)
    str_repr = @inferred mystring(poly, A)
    poly_round = parse_OrePoly(str_repr, A)
    @test poly_round == poly

    mon_str = @inferred mystring(mon(poly, 1), A)
    @test occursin("x", mon_str)

    vec = [poly, parse_OrePoly("dx", A)]
    vec_str = @inferred mystring(vec, A)
    vec_round = parse_vector_OrePoly(vec_str, A)
    @test length(vec_round) == 2
    @test vec_round[1] == poly
    @test vec_round[2] == vec[2]
end
