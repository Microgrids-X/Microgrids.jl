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

function dispatch(Pnl_req,
                  Pbatt_cmax, Pbatt_dmax,
                  Pgen_min , Pgen_max,
                  Pelyz_min, Pelyz_max, 
                  Pfc_min,  Pfc_max)
    Pgen = 0.0
    Pcurt = 0.0
    Pshed = 0.0
    Pelyz = 0.0
    Pfc = 0.0
    # Pnl_req >= 0 - load excess - after evaluating the production (Pnl = Pload - VRE generation)
    if Pnl_req >= 0.0
        Pbatt=min(Pbatt_dmax,Pnl_req)

        if Pnl_req - Pbatt > Pfc_min
            Pfc=min(Pfc_max,Pnl_req-Pbatt)
        end

        if Pnl_req - Pbatt - Pfc > Pgen_min
            Pgen=min(Pgen_max,Pnl_req-Pbatt-Pfc)
        end

        Pshed = Pnl_req - Pbatt - Pfc - Pgen

    # Pnl_req < 0 - VRE excess
    elseif Pnl_req < 0.0
        Pbatt = max(Pbatt_cmax,Pnl_req)
       
        if Pbatt - Pnl_req > Pelyz_min
            Pelyz = min(Pelyz_max,Pbatt-Pnl_req)
        end

        Pcurt = Pbatt - Pelyz - Pnl_req
    end

    Pnl = Pnl_req - Pshed + Pcurt

    return Pnl, Pgen , Pbatt, Pcurt, Pshed, Pelyz, Pfc
end


function dispatch2(Pnl_req,
                    Pbatt_cmax, Pbatt_dmax,
                    Pgen_min , Pgen_max,
                    Pelyz_min, Pelyz_max, 
                    Pfc_min,  Pfc_max)
    Pgen = 0.0
    Pcurt = 0.0
    Pshed = 0.0
    Pelyz = 0.0
    Pfc = 0.0
    # Pnl_req >= 0 - load excess - after evaluating the production (Pnl = Pload - VRE generation)
    if Pnl_req >= 0.0

        if Pnl_req > Pfc_min

            Pfc = min(Pfc_max, Pnl_req)

        end

        Pbatt = min(Pbatt_dmax, Pnl_req - Pfc)

        if Pnl_req - Pbatt - Pfc > Pgen_min
            
            Pgen = min(Pgen_max, Pnl_req - Pbatt - Pfc)

        end

        Pshed = Pnl_req - Pbatt - Pfc - Pgen

    # Pnl_req < 0 - VRE excess
    elseif Pnl_req < 0.0

        if -1 * Pnl_req > Pelyz_min

            Pelyz=min(Pelyz_max, - Pnl_req)

        end
        
        Pbatt=max(Pbatt_cmax,Pnl_req - Pelyz)

        Pcurt = Pbatt - Pelyz - Pnl_req
    end

    Pnl = Pnl_req - Pshed + Pcurt

    return Pnl, Pgen , Pbatt, Pcurt, Pshed, Pelyz, Pfc
end