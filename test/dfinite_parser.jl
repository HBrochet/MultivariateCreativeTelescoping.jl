using Test
using MultivariateCreativeTelescoping

@testset "dfinite_expr_to_ann with provided algebra" begin
    A = OreAlg(order = "lex dx", ratdiffvars = (["x"], ["dx"]))

    gens = dfinite_expr_to_ann("sin(x)", A)
    @test gens[1] == parse_OrePoly("dx^2 + 1", A)

    gens2 = dfinite_expr_to_ann(:(bessel_J(2, x)), A)
    @test gens2[1] == parse_OrePoly("x^2*dx^2 + x*dx + (x^2 - 4)", A)
end

@testset "dfinite_expr_to_ann automatic algebra inference" begin
    gens, A =  dfinite_expr_to_ann("cos(x)")
    @test gens[1] == parse_OrePoly("dx^2 + 1", A)
end

@testset "dfinite_expr_to_ann reports unclosed parenthesis" begin
    err = @test_throws ErrorException dfinite_expr_to_ann("exp(1/(x^2+y^3)")
    @test err.value.msg == "Input expression is incomplete or malformed"
end

@testset "dfinite_expr_to_ann with multiple ratdiff variables" begin
    A = dfinite_ore_alg(["x", "y"])
    gens = dfinite_expr_to_ann("sin(y)", A)
    @test gens[1] == parse_OrePoly("dy^2 + 1", A)
    @test gens[2] == parse_OrePoly("dx", A)
end

@testset "dfinite closure by addition and multiplication" begin
    A = OreAlg(order = "grevlex dx", ratdiffvars = (["x"], ["dx"]))

    gens_add = dfinite_expr_to_ann("sin(x)+cos(x)", A)
    @test gens_add[1] == parse_OrePoly("dx^2 + 1", A)

    gens_sub = dfinite_expr_to_ann("sin(x)-cos(x)", A)
    @test gens_sub[1] == parse_OrePoly("dx^2 + 1", A)

    gens_neg = dfinite_expr_to_ann("-sin(x)", A)
    @test gens_neg[1] == parse_OrePoly("dx^2 + 1", A)

    gens_mul = dfinite_expr_to_ann("sin(x)*cos(x)", A)
    @test gens_mul[1] == parse_OrePoly("dx^3 + 4*dx", A)
end

@testset "dfinite closure by addition and multiplication in 2 variables" begin
    A = dfinite_ore_alg(["x", "y"])
    A = OreAlg(order = "grevlex dx dy", ratdiffvars = (["x","y"], ["dx","dy"]))

    gens_add = dfinite_expr_to_ann("sin(x)+cos(y)", A)
    @test gens_add[1] == parse_OrePoly("dy^3 + dy", A)
    @test gens_add[2] == parse_OrePoly("dx^2 + dy^2 + 1", A)
    @test gens_add[3] == parse_OrePoly("dx*dy", A)

    gens_mul = dfinite_expr_to_ann("sin(x)*cos(y)", A)
    @test length(gens_mul) > 0
    @test gens_mul[1] == parse_OrePoly("dy^2 + 1", A)
    @test gens_mul[2] == parse_OrePoly("dx^2 + 1", A)
end

@testset "dfinite algebraic leaf from minimal polynomial" begin
    A = OreAlg(order = "lex dx", ratdiffvars = (["x"], ["dx"]))
    gens = dfinite_expr_to_ann(:((x + 1)^(1//2)), A)
    @test gens[1] == parse_OrePoly("(2*x+2)*dx - 1", A)
end

@testset "dfinite algebraic leaf in 2 variables" begin
    A = OreAlg(order = "grevlex dx dy", ratdiffvars = (["x", "y"], ["dx", "dy"]))
    gens = dfinite_expr_to_ann(:((x + y + 1)^(1//2)), A)
    @test length(gens) == 2
    @test gens[1] == parse_OrePoly("(2*x + 2*y + 2)*dx + (-1)", A)
    @test gens[2] == parse_OrePoly("(2*x + 2*y + 2)*dy + (-1)", A)
end

@testset "dfinite algebraic leaf with nested square roots" begin
    A = OreAlg(order = "grevlex dx dy", ratdiffvars = (["x", "y"], ["dx", "dy"]))
    gens = dfinite_expr_to_ann(:((x + (y + 1)^(1//2))^(1//2)), A)
    @test length(gens) == 2
    @test gens[1] == parse_OrePoly("2*x*dx + (4*y + 4)*dy - 1", A)
    @test gens[2] == parse_OrePoly("(4*x^2 - 4*y - 4)*dx^2 + 4*x*dx - 1", A)
    @test all(g -> g != parse_OrePoly("dy", A), gens)
end

@testset "dfinite division by hyperexponential denominator" begin
    A = OreAlg(order = "grevlex dx", ratdiffvars = (["x"], ["dx"]))
    gens = dfinite_expr_to_ann(:(sin(x) / exp(x^2)), A)
    @test length(gens) == 1
    @test gens[1] == parse_OrePoly("dx^2 + 4*x*dx + 4*x^2 + 3", A)
end

@testset "dfinite product with composed coefficients" begin
    gens, A = dfinite_expr_to_ann("exp(x^2+y^2+1)*bessel_J(2,x+y)")
    @test gens[1] == parse_OrePoly("(x^2 + 2*x*y + y^2)*dy^2 + (-4*x^2*y - 8*x*y^2 + x - 4*y^3 + y)*dy + (4*x^2*y^2 - x^2 + 8*x*y^3 - 4*x*y + 4*y^4 - 3*y^2 - 4)", A)
    @test gens[2] == parse_OrePoly("(x^2 + 2*x*y + y^2)*dx^2 + (-4*x^3 - 8*x^2*y - 4*x*y^2 + x + y)*dx + (4*x^4 + 8*x^3*y + 4*x^2*y^2 - 3*x^2 - 4*x*y - y^2 - 4)", A)
end


# The annihilator of 1/exp((x^2+y^2+1))*bessel_J(2,x+y) is smaller than ann(bessel_J(2,x+y)/exp((x^2+y^2+1)))
# @testset "dfinite leading constant in hyperexponential reciprocal" begin
#     gens1, A1 = dfinite_expr_to_ann("1/exp((x^2+y^2+1))*bessel_J(2,x+y)")
#     @test A1.inp.ratdiffvars == (["x", "y"], ["dx", "dy"])
#     @test length(gens1) == 2
# end

@testset "dfinite symbolic database parameter in ratvars" begin
    gens1, A1 = dfinite_expr_to_ann("1/exp((x^2+y^2+1))*bessel_J(k,x+y)", ratvars = ["k"])
    A2 = dfinite_ore_alg(["x", "y"]; ratvars = ["k"])
    gens2 = dfinite_expr_to_ann("1/exp((x^2+y^2+1))*bessel_J(k,x+y)", A2)
    @test A1.inp.ratdiffvars == A2.inp.ratdiffvars
    @test A1.inp.ratvars == A2.inp.ratvars
    @test gens1 == gens2
end

@testset "dfinite rational database parameter in ratvars" begin
    gens1, A1 = dfinite_expr_to_ann("bessel_J(1/(k+1),x+y)", ratvars = ["k"])
    @test A1.inp.ratdiffvars == (["x", "y"], ["dx", "dy"])
    @test length(gens1) == 2
end

@testset "dfinite polynomial power with symbolic exponent" begin
    A = dfinite_ore_alg(["x"]; ratvars = ["s"])
    gens = dfinite_expr_to_ann("(x + 1)^s", A)
    @test gens[1] == parse_OrePoly("(x + 1)*dx - s", A)
end

@testset "dfinite polynomial power with symbolic exponent and inferred algebra" begin
    gens1, A1 = dfinite_expr_to_ann("(x + 1)^s", ratvars = ["s"])
    A2 = dfinite_ore_alg(["x"]; ratvars = ["s"])
    gens2 = dfinite_expr_to_ann("(x + 1)^s", A2)
    @test A1.inp.ratvars == ["s"]
    @test A1.inp.ratdiffvars == A2.inp.ratdiffvars
    @test gens1 == gens2
end

@testset "dfinite polynomial power with mixed symbolic exponent reports explicit error" begin
    err = @test_throws ErrorException dfinite_expr_to_ann("x^(s+1/2)", ratvars = ["s"])
    @test err.value.msg == "Polynomial powers currently support only bare symbolic exponents, not mixed expressions like s + 1 / 2"
end

@testset "dfinite product uses polynomial power factor" begin
    A = dfinite_ore_alg(["x"]; ratvars = ["s"])
    gens = dfinite_expr_to_ann("sin(x) * (x + 1)^s", A)
    @test gens[1] == parse_OrePoly("(x^2 + 2*x + 1)*dx^2 + (-2*x*s - 2*s)*dx + (x^2 + 2*x + s^2 + s + 1)", A)
end

#todo: this is bugged @time gens, A  = dfinite_expr_to_ann("sin(x+y)*bessel_J(k,2*x-y)/exp((x^2 + 1)^(1/2))",ratvars=["k"])
