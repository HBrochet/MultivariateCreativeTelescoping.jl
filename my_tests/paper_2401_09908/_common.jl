using MultivariateCreativeTelescoping

function annihilator_from_dfinite_parser(expr::AbstractString, A::OreAlg)
    Ap = to_ore_alg_with_rat_coeffs(A)
    annp = dfinite_expr_to_ann(expr, Ap)
    return map_algebras(annp, Ap, A)
end
