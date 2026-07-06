using Test
using MultivariateCreativeTelescoping
using Nemo

@testset "minimal_polynomial nested radical" begin
    A = OreAlg(order = "grevlex dx", ratdiffvars = (["x"], ["dx"]))
    expr = :((x + (1 + x)^(1//2))^(1//2))

    p = MultivariateCreativeTelescoping.minimal_polynomial(expr, A)
    R = parent(p)
    v = Nemo.gens(R)

    expected = v[end]^4 - 2 * v[1] * v[end]^2 + v[1]^2 - v[1] - 1
    @test iszero(p - expected) || iszero(p + expected)
end

@testset "minimal_polynomial nested radical" begin
    A = OreAlg(order = "grevlex dx", ratdiffvars = (["x"], ["dx"]))
    expr = :(((1 + x)^(1//2) + (1 + x)^(1//2))^(1//2))

    p = MultivariateCreativeTelescoping.minimal_polynomial(expr, A)
    R = parent(p)
    v = Nemo.gens(R)

    expected = v[end]^4 - 4 * v[1] - 4
    @test iszero(p - expected) || iszero(p + expected)
end

@testset "minimal_polynomial two variables" begin
    A = OreAlg(order = "grevlex dx dy", ratdiffvars = (["x", "y"], ["dx", "dy"]))
    expr = :((x + y)^(1//2))

    p = MultivariateCreativeTelescoping.minimal_polynomial(expr, A)
    R = parent(p)
    v = Nemo.gens(R)

    expected = v[end]^2 - v[1] - v[2]
    @test iszero(p - expected) || iszero(p + expected)
end

@testset "minimal_polynomial with division" begin
    A = OreAlg(order = "grevlex dx", ratdiffvars = (["x"], ["dx"]))
    expr = :(((x + 1) / (x + 2))^(1//2))

    p = MultivariateCreativeTelescoping.minimal_polynomial(expr, A)
    R = parent(p)
    v = Nemo.gens(R)

    expected = (v[1] + 2) * v[end]^2 - v[1] - 1
    @test iszero(p - expected) || iszero(p + expected)
end
