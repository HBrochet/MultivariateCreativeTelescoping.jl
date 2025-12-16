using Test
using StaticArrays
using MultivariateCreativeTelescoping

@testset "TypeMonOrder make_order" begin
    dic = Dict("dx" => 1, "dy" => 2, "x" => 3, "y" => 4)
    A = OreMonVE{4, Int16}
    ord = make_order(" grevlex x y dx dy", dic, Val(A))
    @test ord isa AbsMonomialOrder

    m1 = makemon([Int16(1), Int16(0), Int16(0), Int16(0)])
    m2 = makemon([Int16(0), Int16(1), Int16(0), Int16(0)])
    m3 = makemon([Int16(0), Int16(0), Int16(1), Int16(0)])
    m4 = makemon([Int16(0), Int16(0), Int16(0), Int16(1)])

    tab = [m1,m2,m3,m4]
    sort!(tab; order = ord)
    res =  [m1, m2, m4, m3]
    @test tab == [m2, m1, m4, m3]
    @test @inferred lt(ord,m2,m1)
end
