using Test
using Nemo: finite_field, polynomial_ring, fraction_field, numerator, denominator, ZZ, QQ
using Nemo.Native: GF
using MultivariateCreativeTelescoping

@testset "TypeCoefficient: Nmod32Γ arithmetic" begin
    ctx = Nmod32Γ(17)
    a = UInt32(3)
    b = UInt32(5)

    @test @inferred(add(a, b, ctx)) == UInt32(8)
    @test @inferred(sub(a, b, ctx)) == UInt32(15)
    @test @inferred(opp(b, ctx)) == UInt32(12)
    @test @inferred(mul(a, b, ctx)) == UInt32(15)
    @test @inferred(inv(b, ctx)) == UInt32(7)

    prod_buf = @inferred(addmul(UInt64(14), a, b, ctx))
    @test @inferred normal(prod_buf, ctx) == 12
    sub_buf = @inferred(submul(prod_buf, a, b, ctx))
    @test @inferred normal(sub_buf, ctx) == UInt32(14)

    inflated = @inferred(inflate(a, ctx))
    @test inflated == UInt64(a)
    @test @inferred(deflate(inflated, ctx)) == a
    @test @inferred(normal(UInt64(34), ctx)) == UInt32(0)
    @test @inferred(convertn(-3, ctx)) == UInt32(14)
end

@testset "TypeCoefficient: rational univariate mod small p" begin
    p = 7
    Fp = GF(p)
    R, x = polynomial_ring(Fp, "x")
    F = fraction_field(R)
    ctx = UnivRatFunModpCtx(F, R, [x], p)

    f = F(x + 1)
    g = F(x^2 + 1)

    @test @inferred(add(f, g, ctx)) == f + g
    @test @inferred(sub(g, f, ctx)) == g - f
    inv_f = @inferred(inv(f, ctx))
    @test mul(inv_f,f,ctx) == one(F)
    @test @inferred(normal(f, ctx)) === f
    @test @inferred(inflate(f, ctx)) === f
    @test @inferred(deflate(f, ctx)) === f

    @test @inferred(convert(2, ctx)) == ctx.F(2)
    @test @inferred(convertn(-3, ctx)) == ctx.F(-3)

    val = Fp(3)

    @test @inferred(evaluate(f, val)) == Fp(4)
    @test @inferred(evaluate(f, [val])) == Fp(4)
    @test_throws ErrorException evaluate(f, [val, val])
end

@testset "TypeCoefficient: rational univariate mod large p" begin
    p = ZZ(18446744073709551629)
    Fp = GF(p)
    R, x = polynomial_ring(Fp, "x")
    F = fraction_field(R)
    ctx = UnivRatFunModPCtx(F, R, [x], p)

    f = F(x + 1)
    g = F(x^2 + 1)

    @test @inferred(add(f, g, ctx)) == f + g
    @test @inferred(sub(g, f, ctx)) == g - f
    inv_f = @inferred(inv(f, ctx))
    @test @inferred mul(inv_f,f,ctx) == one(F)
    @test @inferred(normal(f, ctx)) === f
    @test @inferred(inflate(f, ctx)) === f
    @test @inferred(deflate(f, ctx)) === f

    @test @inferred(convert(2, ctx)) == ctx.F(2)
    @test @inferred(convertn(-3, ctx)) == ctx.F(-3)

    val = Fp(3)
    @test @inferred(evaluate(f, val)) == F(4)
    @test @inferred(evaluate(f, [val])) == F(4)
    @test_throws ErrorException evaluate(f, [val, val])
end

@testset "TypeCoefficient: rational univariate QQ context" begin
    R, x = polynomial_ring(ZZ, "x")
    F = fraction_field(R)
    ctx = UnivRatFunQQCtx(F, R, [x])

    f = F(x + 1)
    g = F(x^2 + 2)

    @test @inferred(add(f, g, ctx)) == f + g
    @test @inferred(sub(g, f, ctx)) == g - f
    inv_f = @inferred(inv(f, ctx))
    @test @inferred(mul(inv_f,f,ctx)) == one(F)
    @test @inferred(normal(f, ctx)) === f
    @test @inferred(inflate(f, ctx)) === f
    @test @inferred(deflate(f, ctx)) === f

    @test @inferred(convert(5, ctx)) == ctx.F(5)
    @test @inferred(convertn(-4, ctx)) == ctx.F(-4)

    val = ZZ(3)
    @test @inferred(evaluate(f, val)) == ZZ(4)
    @test @inferred(evaluate(f, [val])) == ZZ(4)
    @test_throws ErrorException evaluate(f, [val, val])
end

@testset "TypeCoefficient: rational multivariate mod small p" begin
    p = 11
    Fp = GF(p)
    R, (x, y) = polynomial_ring(Fp, ["x", "y"])
    F = fraction_field(R)
    ctx = RatFunModpCtx(F, R, [x, y], p)

    f = F(x + y)
    g = F(x * y + 1)

    @test @inferred(add(f, g, ctx)) == f + g
    @test @inferred(sub(g, f, ctx)) == g - f
    @test @inferred(mul(f, g, ctx)) == f * g
    @test @inferred(submul(g, f, g, ctx)) == g - f * g
    @test @inferred(opp(f, ctx)) == -f
    inv_f = @inferred(inv(f, ctx))
    @test @inferred(mul(inv_f, f, ctx)) == one(F)
    @test @inferred(normal(f, ctx)) === f
    @test @inferred(inflate(f, ctx)) === f
    @test @inferred(deflate(f, ctx)) === f

    @test @inferred(convert(5, ctx)) == ctx.F(5)
    @test @inferred(convertn(-9, ctx)) == ctx.F(-9)

    vals = [Fp(2), Fp(3)]
    @test @inferred(evaluate(f, vals)) == ctx.F(5)
end

@testset "TypeCoefficient: rational multivariate mod large p" begin
    p = ZZ(18446744073709551629)
    Fp = GF(p)
    R, (x, y) = polynomial_ring(Fp, ["x", "y"])
    F = fraction_field(R)
    ctx = RatFunModPCtx(F, R, [x, y], p)

    f = F(x + y + 1)
    g = F(x * y + 1)

    @test @inferred(add(f, g, ctx)) == f + g
    @test @inferred(sub(g, f, ctx)) == g - f
    @test @inferred(mul(f, g, ctx)) == f * g
    @test @inferred(submul(g, f, g, ctx)) == g - f * g
    @test @inferred(opp(f, ctx)) == -f

    @test @inferred(convert(7, ctx)) == ctx.F(7)
    @test @inferred(convertn(-9, ctx)) == ctx.F(-9)

    vals = [Fp(2), Fp(4)]
    @test @inferred(evaluate(f, vals)) == ctx.F(7)
end

@testset "TypeCoefficient: rational multivariate QQ context" begin
    R, (x, y) = polynomial_ring(ZZ, ["x", "y"])
    F = fraction_field(R)
    ctx = RatFunQQCtx(F, R, [x, y])

    f = F(x + y)
    g = F(x * y + 1)

    @test @inferred(add(f, g, ctx)) == f + g
    @test @inferred(sub(g, f, ctx)) == g - f
    @test @inferred(mul(f, g, ctx)) == f * g
    @test @inferred(submul(g, f, g, ctx)) == g - f * g
    @test @inferred(opp(f, ctx)) == -f

    @test @inferred(convert(5, ctx)) == ctx.F(5)
    @test @inferred(convertn(-9, ctx)) == ctx.F(-9)

    vals = [ZZ(2), ZZ(3)]
    @test @inferred(evaluate(f, vals)) == QQ(5)
end

@testset "TypeCoefficient: QQCtx arithmetic" begin
    ctx = QQCtx()
    a = QQ(2)
    b = QQ(3)

    @test @inferred(add(a, b, ctx)) == QQ(5)
    @test @inferred(sub(b, a, ctx)) == QQ(1)
    @test @inferred(mul(a, b, ctx)) == QQ(6)
    inv_a = @inferred(inv(a, ctx))
    @test inv_a * a == QQ(1)
    @test @inferred(submul(b, a, b, ctx)) == QQ(-3)
    @test @inferred(convert(-4, ctx)) == QQ(-4)
    @test @inferred(convertn(7, ctx)) == QQ(7)
    @test @inferred(normal(a, ctx)) === a
end

# todo I get an error line 198 explain me why
@testset "TypeCoefficient: RingCoeffCtx - ZZRingElem" begin
    ctx = RingCoeffCtx(ZZ)
    a = ZZ(5)
    b = ZZ(-3)
    @test add(a, b, ctx) == ZZ(2)
    @test sub(a, b, ctx) == ZZ(8)
    @test opp(a, ctx) == ZZ(-5)
    @test mul(a, b, ctx) == ZZ(-15)
    @test submul(a, b, a, ctx) == a - b * a
    @test convert(4, ctx) == ZZ(4)
    @test convertn(-6, ctx) == ZZ(-6)
    @test normal(a, ctx) === a
    @test inflate(a, ctx) === a
    @test deflate(a, ctx) === a
    @test zero(ctx) == ZZ(0)
    @test one(ctx) == ZZ(1)
    @test iszero(zero(ctx), ctx)
    @test isone(one(ctx), ctx)
    a1 = deepcopy(a)
    @test add!(a1, b, ctx) == a + b
    b1 = deepcopy(b)
    @test sub!(b1, a, ctx) == b - a
    c1 = deepcopy(a)
    @test mul!(c1, b, ctx) == a * b
    d1 = deepcopy(a)
    @test opp!(d1, ctx) == -a
end

@testset "TypeCoefficient: RingCoeffCtx - ZZPolyRingElem" begin
    R, x = polynomial_ring(ZZ, "x")
    ctx = RingCoeffCtx(R)
    a = x + 1
    b = x^2 - 2x + 3
    @test add(a, b, ctx) == a + b
    @test sub(b, a, ctx) == b - a
    @test opp(b, ctx) == -b
    @test mul(a, b, ctx) == a * b
    @test submul(b, a, b, ctx) == b - a * b
    @test convert(7, ctx) == R(7)
    @test convertn(-4, ctx) == R(-4)
    @test normal(a, ctx) === a
    @test inflate(b, ctx) === b
    @test deflate(a, ctx) === a
    @test zero(ctx) == R(0)
    @test one(ctx) == R(1)
    @test iszero(zero(ctx), ctx)
    @test isone(one(ctx), ctx)
    a1 = deepcopy(a)
    @test add!(a1, b, ctx) == a + b
    b1 = deepcopy(b)
    @test sub!(b1, a, ctx) == b - a
    c1 = deepcopy(a)
    @test mul!(c1, b, ctx) == a * b
    d1 = deepcopy(a)
    @test opp!(d1, ctx) == -a
end

@testset "TypeCoefficient: RingCoeffCtx - ZZMPolyRingElem" begin
    R, (x, y) = polynomial_ring(ZZ, ["x", "y"])
    ctx = RingCoeffCtx(R)
    a = x + y
    b = x * y + 2
    @test add(a, b, ctx) == a + b
    @test sub(b, a, ctx) == b - a
    @test opp(a, ctx) == -a
    @test mul(a, b, ctx) == a * b
    @test submul(a, b, a, ctx) == a - b * a
    @test convert(3, ctx) == R(3)
    @test convertn(-5, ctx) == R(-5)
    @test normal(a, ctx) === a
    @test inflate(b, ctx) === b
    @test deflate(a, ctx) === a
    @test zero(ctx) == R(0)
    @test one(ctx) == R(1)
    a1 = deepcopy(a)
    @test add!(a1, b, ctx) == a + b
    b1 = deepcopy(b)
    @test sub!(b1, a, ctx) == b - a
    c1 = deepcopy(a)
    @test mul!(c1, b, ctx) == a * b
    d1 = deepcopy(a)
    @test opp!(d1, ctx) == -a
end

@testset "TypeCoefficient: RingCoeffCtx - fpPolyRingElem" begin
    R, x = polynomial_ring(GF(7), "x")
    ctx = RingCoeffCtx(R)
    a = x + 2
    b = x^2 + 3
    @test add(a, b, ctx) == a + b
    @test sub(b, a, ctx) == b - a
    @test opp(a, ctx) == -a
    @test mul(a, b, ctx) == a * b
    @test submul(b, a, b, ctx) == b - a * b
    @test convert(5, ctx) == R(5)
    @test convertn(-3, ctx) == R(-3)
    @test normal(a, ctx) === a
    @test inflate(b, ctx) === b
    @test deflate(a, ctx) === a
    @test zero(ctx) == R(0)
    @test one(ctx) == R(1)
    a1 = deepcopy(a)
    @test add!(a1, b, ctx) == a + b
    b1 = deepcopy(b)
    @test sub!(b1, a, ctx) == b - a
    c1 = deepcopy(a)
    @test mul!(c1, b, ctx) == a * b
    d1 = deepcopy(a)
    @test opp!(d1, ctx) == -a
end

@testset "TypeCoefficient: RingCoeffCtx - fpMPolyRingElem" begin
    R, (x, y) = polynomial_ring(GF(5), ["x", "y"])
    ctx = RingCoeffCtx(R)
    a = x + y + 1
    b = x * y + 2
    @test add(a, b, ctx) == a + b
    @test sub(b, a, ctx) == b - a
    @test opp(b, ctx) == -b
    @test mul(a, b, ctx) == a * b
    @test submul(a, b, a, ctx) == a - b * a
    @test convert(4, ctx) == R(4)
    @test convertn(-2, ctx) == R(-2)
    @test normal(a, ctx) === a
    @test inflate(b, ctx) === b
    @test deflate(a, ctx) === a
    a1 = deepcopy(a)
    @test add!(a1, b, ctx) == a + b
    b1 = deepcopy(b)
    @test sub!(b1, a, ctx) == b - a
    c1 = deepcopy(a)
    @test mul!(c1, b, ctx) == a * b
    d1 = deepcopy(a)
    @test opp!(d1, ctx) == -a
end

@testset "TypeCoefficient: RingCoeffCtx - FpPolyRingElem" begin
    p = ZZ(2305843009213693951) # large prime to force Fp type
    Fp = GF(p)
    R, x = polynomial_ring(Fp, "x")
    ctx = RingCoeffCtx(R)
    a = x + 1
    b = x^2 + Fp(3)
    @test add(a, b, ctx) == a + b
    @test sub(b, a, ctx) == b - a
    @test opp(a, ctx) == -a
    @test mul(a, b, ctx) == a * b
    @test submul(b, a, b, ctx) == b - a * b
    @test convert(6, ctx) == R(6)
    @test convertn(-7, ctx) == R(-7)
    @test normal(a, ctx) === a
    @test inflate(b, ctx) === b
    @test deflate(a, ctx) === a
    a1 = deepcopy(a)
    @test add!(a1, b, ctx) == a + b
    b1 = deepcopy(b)
    @test sub!(b1, a, ctx) == b - a
    c1 = deepcopy(a)
    @test mul!(c1, b, ctx) == a * b
    d1 = deepcopy(a)
    @test opp!(d1, ctx) == -a
end

@testset "TypeCoefficient: RingCoeffCtx - FpMPolyRingElem" begin
    p = ZZ(2305843009213693951)
    Fp = GF(p)
    R, (x, y) = polynomial_ring(Fp, ["x", "y"])
    ctx = RingCoeffCtx(R)
    a = x + y + 1
    b = x * y + Fp(2)
    @test add(a, b, ctx) == a + b
    @test sub(b, a, ctx) == b - a
    @test opp(b, ctx) == -b
    @test mul(a, b, ctx) == a * b
    @test submul(a, b, a, ctx) == a - b * a
    @test convert(9, ctx) == R(9)
    @test convertn(-11, ctx) == R(-11)
    @test normal(a, ctx) === a
    @test inflate(b, ctx) === b
    @test deflate(a, ctx) === a
    a1 = deepcopy(a)
    @test add!(a1, b, ctx) == a + b
    b1 = deepcopy(b)
    @test sub!(b1, a, ctx) == b - a
    c1 = deepcopy(a)
    @test mul!(c1, b, ctx) == a * b
    d1 = deepcopy(a)
    @test opp!(d1, ctx) == -a
end
