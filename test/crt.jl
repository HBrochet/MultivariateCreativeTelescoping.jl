@testset "CRT With/Without Tracer (F4)" begin

    Rpt = OreAlg(order = "grevlex p1 p2 p3 p4 > grevlex d1 d2 d3 d4",
                 poldiffvars = (["p1", "p2", "p3", "p4"], ["d1", "d2", "d3", "d4"]))
    Gstring = "[(-7)*p1^3*1 + (21)*d1*p1^2*1 + (21)*p1*p2*1 + (-21)*d1^2*p1*1 + (42)*d2*p1*1 + (1)*p1*1 + (-21)*d1*p2*1 + (-14)*p3*1 + (7)*d1^3*1 + (-42)*d1*d2*1 + (42)*d3*1, (-21)*p1^2*1 + (42)*d1*p1*1 + (22)*p2*1 + (-21)*d1^2*1 + (42)*d2*1, (-42)*p1*1 + (1)*p3*1 + (42)*d1*1, (1)*p4*1 + (-42)*1]"
    G = parse_vector_OrePoly(Gstring, Rpt)

    res = compute_with_CRT(f4, Rpt, G, Rpt)

    s= "[(1804)*p2^2*1 + (7140)*d2*p2*1 + (-73914)*p2*1 + (3)*d1*p3*1 + (126)*d3*p3*1 + (7056)*d2^2*1 + (-140868)*d2*1 + (-75768)*1, (41)*p2*p3*1 + (84)*d2*p3*1 + (-1719)*p3*1 + (126)*d1*1 + (5292)*d3*1, (1)*p3^2*1 + (-1848)*p2*1 + (-3528)*d2*1 + (-1764)*1, (42)*p1*1 + (-1)*p3*1 + (-42)*d1*1, (1)*p4*1 + (-42)*1]"
    res0 = parse_vector_OrePoly(s,Rpt)
    @test res == res0 

    par_tr = f4_param(tracer = Val(:learn))
    res_tr = compute_with_CRT(f4, Rpt, G, Rpt, par_tr, param = crt_param(tracer = Val(true)))

    @test res_tr == res0
end
