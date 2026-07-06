using Test
using MultivariateCreativeTelescoping

@testset "ann_poly_power direct one variable" begin
    A = dfinite_ore_alg(["x"]; ratvars = ["s"])
    gens = dfinite_expr_to_ann(:((x + 1)^s), A)
    @test length(gens) == 1
    @test gens[1] == parse_OrePoly("(x + 1)*dx - s", A)
end

@testset "ann_poly_power direct two variables" begin
    A = dfinite_ore_alg(["x", "y"]; ratvars = ["s"])
    gens = dfinite_expr_to_ann("((x + y + 1)^s)", A)
    @test length(gens) == 2
    @test parse_OrePoly("(x + y + 1)*dx - s", A) in gens
    @test parse_OrePoly("(x + y + 1)*dy - s", A) in gens
end
