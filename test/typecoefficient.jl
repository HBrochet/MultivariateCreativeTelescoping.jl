using Test
using Nemo: finite_field, polynomial_ring, fraction_field, numerator, denominator, ZZ, QQ, divexact
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

	@testset "TypeCoefficient: ring contexts (fraction-free)" begin
	    function test_ringctx_ops(
	        ctx::RingCtx{T,Tbuf},
	        a::T,
	        b::T,
	        c::T,
	        div_a::T,
	        div_b::T,
	        gcd_a::T,
	        gcd_b::T,
	    ) where {T,Tbuf}
	        @test @inferred(opp(a, ctx)) == -a
	        @test @inferred(add(a, b, ctx)) == a + b
	        @test @inferred(sub(a, b, ctx)) == a - b
	        @test @inferred(mul(a, b, ctx)) == a * b
	
	        aa = deepcopy(a)
	        @test @inferred(opp!(aa, ctx)) == -a
	
	        aa = deepcopy(a)
	        @test @inferred(add!(aa, b, ctx)) == a + b
	
	        aa = deepcopy(a)
	        @test @inferred(sub!(aa, b, ctx)) == a - b
	
	        aa = deepcopy(a)
	        @test @inferred(mul!(aa, b, ctx)) == a * b
	
	        aa = deepcopy(div_a)
	        @test @inferred(divexact!(aa, div_b, ctx)) == divexact(div_a, div_b)

	        @test @inferred(submul(a, b, c, ctx)) == a - b * c

        inflated = @inferred(inflate(a, ctx))
        @test @inferred(normal(inflated, ctx)) === inflated
        @test @inferred(deflate(inflated, ctx)) == a

        @test @inferred(Base.gcd(gcd_a, gcd_b, ctx)) == Base.gcd(gcd_a, gcd_b)

        @test @inferred(zero(ctx)) == ctx.R(0)
        @test @inferred(one(ctx)) == ctx.R(1)
        @test @inferred(zero(T, ctx)) == ctx.R(0)
        @test @inferred(one(T, ctx)) == ctx.R(1)
        @test @inferred(iszero(zero(ctx), ctx))
        @test @inferred(isone(one(ctx), ctx))

        @test @inferred(convert(5, ctx)) == ctx.R(5)
        @test @inferred(convertn(-3, ctx)) == ctx.R(-3)

        @test_throws ErrorException inv(a, ctx)
    end

	    @testset "ZZCtx" begin
	        ctx = zz_ctx()
	        a = ZZ(6)
	        b = ZZ(4)
	        c = ZZ(3)
	        div_a = ZZ(12)
	        div_b = ZZ(3)
	        gcd_a = ZZ(12)
	        gcd_b = ZZ(18)
	        test_ringctx_ops(ctx, a, b, c, div_a, div_b, gcd_a, gcd_b)
	    end

	    @testset "ZZPolyCtx" begin
	        ctx = zzpoly_ctx("x")
	        x = ctx.vars[1]
	        a = x^2 + 2*x + 1
	        b = x + 3
	        c = x + 1
	        div_a = (x + 1) * (x + 2)
	        div_b = x + 1
	        gcd_a = (x + 1) * (x + 2)
	        gcd_b = (x + 1) * (x + 5)
	        test_ringctx_ops(ctx, a, b, c, div_a, div_b, gcd_a, gcd_b)
	    end

	    @testset "ZZMPolyCtx" begin
	        ctx = zzmpoly_ctx(["x", "y"])
	        x, y = ctx.vars
	        a = x^2 + y + 1
	        b = x*y + 2
	        c = x + y
	        div_a = (x + y) * (x + 1)
	        div_b = x + y
	        gcd_a = (x + y) * (x + 1)
	        gcd_b = (x + y) * (y + 2)
	        test_ringctx_ops(ctx, a, b, c, div_a, div_b, gcd_a, gcd_b)
	    end

	    @testset "fpPolyCtx" begin
	        ctx = fppoly_ctx(7, "x")
	        x = ctx.vars[1]
	        a = x^3 + 2*x + 1
	        b = x + 3
	        c = x + 1
	        div_a = (x + 1) * (x + 2)
	        div_b = x + 1
	        gcd_a = (x + 1) * (x + 2)
	        gcd_b = (x + 1) * (x + 5)
	        test_ringctx_ops(ctx, a, b, c, div_a, div_b, gcd_a, gcd_b)
	    end

	    @testset "fpMPolyCtx" begin
	        ctx = fpmpoly_ctx(7, ["x", "y"])
	        x, y = ctx.vars
	        a = x^2 + y^2 + 1
	        b = x + y
	        c = x + 2*y + 1
	        div_a = (x + y) * (x + 1)
	        div_b = x + y
	        gcd_a = (x + y) * (x + 1)
	        gcd_b = (x + y) * (y + 2)
	        test_ringctx_ops(ctx, a, b, c, div_a, div_b, gcd_a, gcd_b)
	    end
	end
