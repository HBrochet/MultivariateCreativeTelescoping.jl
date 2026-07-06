using Test
using MultivariateCreativeTelescoping

@testset "module Groebner basis from dfinite annihilator" begin
    ann, A0 = dfinite_expr_to_ann("(x + 1)^s", ratvars = ["s"])
    res = compute_module_groebner_basis(ann, A0, ["s"], ["x"], 2)

    @test res.module_var_names == ["ed0", "edx", "edx2"]
    @test length(res.linear_generators) == 2
    @test length(res.square_zero_generators) == 6
    @test !isempty(res.basis)

    x = res.variables["x"]
    s = res.variables["s"]
    ed0 = res.variables["ed0"]
    edx = res.variables["edx"]
    edx2 = res.variables["edx2"]

    @test res.linear_generators[1] == (x + 1) * edx - s * ed0
    @test res.linear_generators[2] == (x + 1) * edx2 + (1 - s) * edx
end

@testset "heuristic Weyl closure via module Groebner basis" begin
    ann, A0 = dfinite_expr_to_ann("(x + 1)^s", ratvars = ["s"])
    res = heuristic_weyl_closure_holonomic_ideal(ann, A0)

    @test res.holonomic
    @test res.trunc_order == 1
    @test res.module_data.localization_var_name == "T"
    @test isholonomic(res.ideal, res.algebra)

    expected = parse_OrePoly("(x + 1)*dx - s", res.algebra)
    @test any(==(expected), res.generators)
end

@testset "heuristic Weyl closure without user parameters" begin
    ann, A0 = dfinite_expr_to_ann("sin(x)")
    res = heuristic_weyl_closure_holonomic_ideal(ann, A0; localizing_factor = :none)

    @test res.holonomic
    @test res.trunc_order == 2
    @test isholonomic(res.ideal, res.algebra)

    expected = parse_OrePoly("dx^2 + 1", res.algebra)
    @test any(==(expected), res.generators)
end
