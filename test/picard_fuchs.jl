using Test
using MultivariateCreativeTelescoping

@testset "picard_fuchs type inference without debug" begin
    MCT = MultivariateCreativeTelescoping
    num, den, A, finalA = MCT._pf_parse_and_homogenize("1/(1-t*x)", ["t"], String[])

    result = @inferred MCT._picard_fuchs(num, den, A, finalA; rho = 1)

    expected = parse_OrePoly("t*dt + 1", result[2])
    @test result[1][1] == expected
end
