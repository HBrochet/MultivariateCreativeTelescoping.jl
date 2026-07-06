using JET
using MultivariateCreativeTelescoping
using Cthulhu 

Rpt = OreAlg(order = "grevlex p1 p2 p3 p4 > grevlex d1 d2 d3 d4",
                poldiffvars = (["p1", "p2", "p3", "p4"], ["d1", "d2", "d3", "d4"]))
Gstring = "[(-7)*p1^3*1 + (21)*d1*p1^2*1 + (21)*p1*p2*1 + (-21)*d1^2*p1*1 + (42)*d2*p1*1 + (1)*p1*1 + (-21)*d1*p2*1 + (-14)*p3*1 + (7)*d1^3*1 + (-42)*d1*d2*1 + (42)*d3*1, (-21)*p1^2*1 + (42)*d1*p1*1 + (22)*p2*1 + (-21)*d1^2*1 + (42)*d2*1, (-42)*p1*1 + (1)*p3*1 + (42)*d1*1, (1)*p4*1 + (-42)*1]"
G = parse_vector_OrePoly(Gstring, Rpt)
par = f4_param()

@code_warntype f4(G, Rpt; param = par)
@report_opt f4(G, Rpt; param = par)

par = F4Param()
@descend f4(G, Rpt; param = par)
@descend_code_warntype f4(G, Rpt)

# todo: F4Param is type unstable ? 