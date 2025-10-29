using Test
using StaticArrays
using MultivariateCreativeTelescoping

@testset "TypeOreMonomial basics" begin
    m1 = makemon([1, 0, 2])
    m2 = makemon([0, 3, 1])

    @test exp(m1) == SVector(1, 0, 2)
    @test getexp(m1, 2) == 0
    @test degree(m1) == 3
    @test length(m1) == 3
    @test nvars(m1) == 3
    @test nvars(typeof(m1)) == 3
    @test exptype(typeof(m1)) == Int
    @test m1[1] == 1
    @test sum(m1) == 3

    prod = @inferred m1 * m2
    @test prod == makemon([1, 3, 3])

    quot = @inferred prod / m1
    @test quot == m2

    pow = @inferred m1^2
    @test pow == makemon([2, 0, 4])

    l = @inferred lcm(m1, m2)
    @test l == makemon([1, 3, 2])
end
