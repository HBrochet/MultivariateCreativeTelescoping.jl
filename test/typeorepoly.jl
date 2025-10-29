using Test
using MultivariateCreativeTelescoping

@testset "TypeOrePolynomial structure" begin
    m1 = makemon([1, 0])
    m2 = makemon([0, 1])
    p = OrePoly([2, 3], [m1, m2])

    @test mons(p) == [m1, m2]
    @test mon(p, 1) === m1
    @test lm(p) === m1
    @test coeffs(p) == [2, 3]
    @test coeff(p, 2) == 3
    @test lc(p) == 2

    @test length(p) == 2
    @test size(p) == (2,)
    q = copy(p)
    @test q !== p
    @test coeffs(q) == coeffs(p)
    @test !iszero(p)
    z = zero(p)
    @test iszero(z)

    @test p[1] == (2, m1)

    resize!(p, 1)
    @test length(p) == 1
    @test coeff(p, 1) == 2

    p[1] = (5, m2)
    @test coeff(p, 1) == 5
    @test mon(p, 1) === m2

    resize!(p, 2)
    p[2] = (7, m1)
    @test degree(p) == 1
end

@testset "ReuseOrePoly behaviour" begin
    buffer_poly = OrePoly([0, 0], [makemon([0, 0]), makemon([0, 0])])
    rp = ReuseOrePoly(buffer_poly, 0)

    @test length(rp) == 2
    @test iszero(rp)

    push!(rp, 3, makemon([1, 0]))
    push!(rp, 4, makemon([0, 1]))
    @test rp.ind == 2
    @test rp[1] == (3, makemon([1, 0]))
    @test rp[2] == (4, makemon([0, 1]))

    push!(rp, 5, makemon([2, 0])) # triggers grow_to!
    @test rp.ind == 3
    @test rp[3] == (5, makemon([2, 0]))
    @test length(rp.op.coeffs) == 4

    rp[2] = (6, makemon([1, 1]))
    @test rp[2] == (6, makemon([1, 1]))
end
