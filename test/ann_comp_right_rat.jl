using Test
using MultivariateCreativeTelescoping

@testset "ann_comp_right_rat univariate closure" begin
    A0 = OreAlg(order = "grevlex dx", ratdiffvars = (["x"], ["dx"]))
    A = _dfinite_temp_closure_alg(A0)

    g = :(x^2 + 1)
    Lt = ann_comp_right_rat(:exp, g, A)[1]

    @test Lt == parse_OrePoly("1/(2*x)*dx - 1", A)
end

@testset "ann_comp_right_rat ignores second closure block" begin
    A0 = dfinite_ore_alg(["x", "y"])
    A = _dfinite_temp_closure_alg(A0)

    g = :(x + y)
    Lt = ann_comp_right_rat(:exp, g, A)

    @test Lt == parse_vector_OrePoly("[dx - 1, dy - 1]", A)
end

@testset "dfinite_expr_to_ann ratfun composition branch" begin
    A = OreAlg(order = "grevlex dx", ratdiffvars = (["x"], ["dx"]))
    gens = dfinite_expr_to_ann(:(exp(x^2 + 1)), A)

    @test gens == [parse_OrePoly("1/(2*x)*dx - 1", A)]
end
