using Test
using MultivariateCreativeTelescoping

@testset "picard_fuchs type inference without debug" begin
    MCT = MultivariateCreativeTelescoping
    num, den, A, finalA = MCT._pf_parse_and_homogenize("1/(1-t*x)", ["t"], String[])

    result = @inferred MCT._picard_fuchs(num, den, A, finalA; rho = 1)

    expected = parse_OrePoly("t*dt + 1", result[2])
    @test result[1][1] == expected
end

@testset "PicardFuchs Rham-Koszul reduction" begin
    data = PicardFuchs("1/(1-t*x)", variables = ["x"], return_data = true, max_order = 4)

    expected = parse_OrePoly("t*dt + 1", data.operator_algebra)
    @test data.operator == expected
    @test all(i in data.homogeneous_algebra.nomul for i in data.homogeneous_algebra.nrdv+1:data.homogeneous_algebra.nrdv+data.homogeneous_algebra.npdv)
    @test !isempty(data.groebner_basis)

    with_param = PicardFuchs("1/(m-t*x)", variables = ["x"], ratvars = ["m"], return_data = true, max_order = 4)
    @test with_param.operator == parse_OrePoly("t*dt + 1", with_param.operator_algebra)
    @test haskey(with_param.homogeneous_algebra.ratvars, "m")
end

@testset "PicardFuchs Alin examples" begin
    s2 = "1/((v2 - 1)*(v1 + v2)*t - v1*v2*(v1 + v2 - 1))"
    data2 = PicardFuchs(s2; variables = ["v1", "v2"], return_data = true, max_order = 2, max_reduction_order = 2)
    expected2 = parse_OrePoly(
        "(t^3 + 11*t^2 - t)*dt^2 + (3*t^2 + 22*t - 1)*dt + (t + 3)",
        data2.operator_algebra,
    )
    @test data2.operator == expected2

    s3 = "1/((v2 + v3 - 1)*(v1 + v3)*(v1 + v2)*t - v3*v2*v1*(v1 + v2 + v3 - 1))"
    data3 = PicardFuchs(s3; variables = ["v1", "v2", "v3"], return_data = true, max_order = 3, max_reduction_order = 2)
    expected3 = parse_OrePoly(
        "(16*t^4 - 24*t^3 + t^2)*dt^3 + (96*t^3 - 108*t^2 + 3*t)*dt^2 + (112*t^2 - 80*t + 1)*dt + (16*t - 4)",
        data3.operator_algebra,
    )
    @test data3.operator == expected3
    @test data3.reduction_order == 2
end
