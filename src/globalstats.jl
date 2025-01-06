struct Statistics
    timings :: Dict{Symbol, Float64}
    counters :: Dict{Symbol, Int}
end

const globalstats = Statistics(Dict(), Dict())

function initglobalstats!()
    for p in []
        push!(globalstats.timings, p => 0.0)
    end

    for p in [:f4_eliminated_spairs_with_prod_crit, :f4_candidate_spairs, :f4_field_operations,:f4_line_reductions,:f4_eliminated_spairs_with_GM,
              :f5_eliminated_signatures, :f5_candidate_signatures, :f5_field_operations,:f5_number_divisions,:f5_size_coeff,:number_primes,:number_evaluation,
              :f5_eliminated_signatures_stophol,:number_pol_evaluation,:f5_size_m,:f5_m_is_pol]
        push!(globalstats.counters, p => 0)
    end
end

initglobalstats!()


