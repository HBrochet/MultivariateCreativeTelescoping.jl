struct Statistics
    timings :: Dict{Symbol, Float64}
    counters :: Dict{Symbol, Int}
    reducers :: Dict{Symbol,Int}
end

const globalstats = Statistics(Dict(), Dict(),Dict())

function initglobalstats!()
    for p in []
        push!(globalstats.timings, p => 0.0)
    end

    for p in [:bbg_nb_spair, :bbg_nb_div, :bbg_size_reducer,:bbg_size_m,:bbg_deg_reducer,
              :f4_eliminated_spairs_with_prod_crit, :f4_candidate_spairs, :f4_field_operations,:f4_line_reductions,:f4_eliminated_spairs_with_GM,
              :f4_nb_reducer_computed, :f4_nb_reducer_used,:f4_size_reducer,:f4_size_m,:f4_deg_reducer,
              :f5_eliminated_signatures, :f5_candidate_signatures, :f5_size_reducer,:f5_deg_reducer,:f5_number_divisions,
              :f5_eliminated_signatures_stophol,:f5_size_m,
              :number_primes,:number_evaluation,:number_pol_evaluation,
              :mct_red_zero,
              :mri_nb_points,
              :f4_max_col_size, :f4_max_row_size
              ]
        push!(globalstats.counters, p => 0)
    end
    empty!(globalstats.reducers)
end

initglobalstats!()


