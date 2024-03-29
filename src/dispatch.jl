"""
    dispatch(Pnl_req, Pbatt_cmax, Pbatt_dmax, Pgen_max)

Return the dispatch decision for a rule-based strategy.

# Arguments
- `Pnl_req`: the requested load at time t.
- `Pbatt_cmax`: the maximum battery charge power at time t.
- `Pbatt_dmax`: the maximum battery discharge power at time t.
- `Pgen_max`: the diesel generator rated power.
"""
function dispatch(Pnl_req,
                  Pbatt_cmax,
                  Pbatt_dmax,
                  Pgen_max)
    Pcurt = 0.0
    Pshed = 0.0
    # Pnl_req >= 0 - load excess - after evaluating the production (Pnl = Pload - VRE generation)
    if Pnl_req >= 0.0
        # battery discharging --> Pbatt positive
        if Pnl_req >= Pbatt_dmax    # max(battery)
            Pbatt = Pbatt_dmax      # max(battery)
            # TODO - also test if (Pnl_req - Pbatt) < dieselgenerator.minimum_load_ratio * gen.power_rated --> usar igual a zero
            # can we use the excedent to charge the battery?
            # dispatch battery or generator based in the fixed+marginal cost?
            if Pnl_req - Pbatt >= Pgen_max  # max(generator)
                Pgen  = Pgen_max
                Pshed = Pnl_req - Pbatt - Pgen
            else
                Pgen = Pnl_req - Pbatt
            end
        else
            Pgen  = 0.0
            Pbatt = Pnl_req
        end
    # Pnl_req < 0 - VRE excess
    elseif Pnl_req < 0.0
        Pgen  = 0.0
        # battery charging --> Pbatt negative
        if Pnl_req >= Pbatt_cmax    # min(battery)
            Pbatt = Pnl_req
        else
            Pbatt = Pbatt_cmax      # min(battery)
            Pcurt  = Pbatt - Pnl_req
        end
    end
    Pnl = Pnl_req - Pshed + Pcurt
    return Pnl, Pgen, Pbatt, Pcurt, Pshed
end

