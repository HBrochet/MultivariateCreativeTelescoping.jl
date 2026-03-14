using Test
using MultivariateCreativeTelescoping

@testset "ann_poly_power direct one variable" begin
    A = dfinite_ore_alg(["x"]; ratvars = ["s"])
    gens = dfinite_expr_to_ann(:((x + 1)^s), A)
    @test length(gens) == 1
    @test gens[1] == parse_OrePoly("dx - s//(x + 1)", A)
end

@testset "ann_poly_power direct two variables" begin
    A = dfinite_ore_alg(["x", "y"]; ratvars = ["s"])
    gens = dfinite_expr_to_ann("((x + y + 1)^s)", A)
    @test length(gens) == 2
    @test parse_OrePoly("dx - s//(x + y + 1)", A) in gens
    @test parse_OrePoly("dy - s//(x + y + 1)", A) in gens
end