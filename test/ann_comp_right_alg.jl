using Test
using MultivariateCreativeTelescoping

@testset "dfinite algebraic composition" begin
    A = OreAlg(order = "grevlex dx dy dz", ratdiffvars = (["x", "y", "z"], ["dx", "dy", "dz"]))
    gens = dfinite_expr_to_ann(:(exp((x + y + 1)^(1//2))), A)

    @test parse_OrePoly("-dx + dy", A) in gens
    @test parse_OrePoly("dx^2 + 1//(2*x + 2*y + 2)*dx - 1//(4*x + 4*y + 4)", A) in gens
    @test parse_OrePoly("dz", A) in gens
end
