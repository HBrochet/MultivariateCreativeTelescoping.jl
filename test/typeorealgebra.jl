using Test
using Nemo: gens
using MultivariateCreativeTelescoping

@testset "TypeOreAlgebra basics" begin
    A = OreAlg(order = "lex x dx", poldiffvars = (["x"], ["dx"]))

    idx_x = A.strvar_to_indexp["x"]
    idx_dx = A.strvar_to_indexp["dx"]

    mx = @inferred makemon(idx_x, A)
    mdx = @inferred makemon(idx_dx, A)

    c1 = @inferred convert(1, ctx(A))
    c2 = @inferred convert(2, ctx(A))

    px = @inferred makepoly(c1, mx)
    pdx = @inferred makepoly(c2, mdx)

    @test maxdeg(px) == 1
    @test maxdeg([px, pdx], idx_x) == 1
    @test @inferred maxdeg(px, idx_x) == 1

    @test @inferred divide(mx, mx)
    @test !divide(mx, mdx)
    @test @inferred divide(mx, mx, A)
    @test !divide(mx, mdx, A)

    @test @inferred iscompatible(mx, mx, A)

    zero_poly = @inferred zero(A)
    one_poly = @inferred one(A)
    @test iszero(zero_poly)
    @test !iszero(one_poly)

    undef = @inferred undefOrePoly(2, A)
    @test length(undef.coeffs) == 2
    @test length(undef.mons) == 2
end

@testset "TypeOreAlgebra fraction-free coefficients" begin
    A_frac = OreAlg(order = "lex x dx", ratvars = ["t"], poldiffvars = (["x"], ["dx"]))
    @test ctx(A_frac) isa UnivRatFunQQCtx

    A_poly = OreAlg(order = "lex x dx", ratvars = ["t"], poldiffvars = (["x"], ["dx"]), fraction_free = true)
    @test ctx(A_poly) isa RingCoeffCtx
    @test A_poly.ratvars["t"] == first(collect(gens(ctx(A_poly).R)))

     A_poly = OreAlg(order = "lex x dt dx", ratdiffvars = (["t"],["dt"]), poldiffvars = (["x"], ["dx"]), fraction_free = true)
    @test ctx(A_poly) isa RingCoeffCtx
end
