"""
    dispatch_1(Pnl_req,
                  Pbatt_cmax,
                  Pbatt_dmax,
                  Pgen_max)

Return the dispatch decision for a rule-based strategy.

# Arguments
- `Pnl_req`: the requested load at time t.
- `Pbatt_cmax`: the maximum battery charge power at time t.
- `Pbatt_dmax`: the maximum battery discharge power at time t.
- `Pgen_max`: the diesel generator rated power.
"""
function dispatch_1(Pnl_req,
                  Pbatt_cmax,
                  Pbatt_dmax,
                  Pgen_max)
    Pspill = 0.0
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
            Pspill  = Pbatt - Pnl_req
        end
    end
    Pnl = Pnl_req - Pshed + Pspill
    return Pnl, Pgen, Pbatt, Pspill, Pshed
end

"""
    dispatch_1(Pnl_req,
                  Pbatt_cmax, Pbatt_dmax,
                  Pgen_min , Pgen_max,
                  Pelyz_min, Pelyz_max, 
                  Pfc_min,  Pfc_max)

# Arguments
- `Pnl_req`: the requested load at time t.
- `Pbatt_cmax`: the maximum battery charge power at time t.
- `Pbatt_dmax`: the maximum battery discharge power at time t.
- `Pgen_max`: the diesel generator rated power.
- `Pelyz_min` : the minimum electrolyzer power at time t.
- `Pelyz_max` : the maximum electrolyzer power at time t.
- 'Pfc_min' : the minimum electrolyzer power at time t.
- 'Pfc_max' : the minimum electrolyzer power at time t.
"""
function dispatch_1(Pnl_req,
                  Pbatt_cmax, Pbatt_dmax,
                  Pgen_min , Pgen_max,
                  Pelyz_min, Pelyz_max, 
                  Pfc_min,  Pfc_max)
    Pgen = 0.0
    Pspill = 0.0
    Pshed = 0.0
    Pelyz = 0.0
    Pfc = 0.0
    
    if Pnl_req >= 0.0
        ## Pnl_req >= 0 - load excess - after evaluating the production (Pnl = Pload - VRE generation)
        # battery discharging --> Pbatt positive
        Pbatt = min(Pbatt_dmax,Pnl_req)

        if Pnl_req - Pbatt > Pfc_min
            # fuel_cell starts on if the battery can't supply the power net load request
            Pfc = min(Pfc_max,Pnl_req-Pbatt)
        end
        
        if Pnl_req - Pbatt - Pfc > Pgen_min
            # Diesel generator starts on if the battery & the fuel_cell can't supply the power net load request
            Pgen = min(Pgen_max,Pnl_req-Pbatt-Pfc)
        end

        Pshed = Pnl_req - Pbatt - Pfc - Pgen

    ## Pnl_req < 0 - VRE excess
    elseif Pnl_req < 0.0
        # battery charging --> Pbatt negative
        Pbatt = max(Pbatt_cmax,Pnl_req)
       
        if Pbatt - Pnl_req > Pelyz_min
            #Electrolyzer starts on if its remain enough power    
            Pelyz = min(Pelyz_max,Pbatt-Pnl_req)
        end

        Pspill = Pbatt - Pelyz - Pnl_req
    end

    Pnl = Pnl_req - Pshed + Pspill

    return Pnl, Pgen , Pbatt, Pspill, Pshed, Pelyz, Pfc
end

"""
    dispatch_2(Pnl_req,
                    Pbatt_cmax, Pbatt_dmax,
                    Pgen_min , Pgen_max,
                    Pelyz_min, Pelyz_max, 
                    Pfc_min,  Pfc_max)
This dispatch strategy use H2 chain components (Fuel_cell and Electrolyzer) prior, and after, the battery.

# Arguments
- `Pnl_req`: the requested load at time t.
- `Pbatt_cmax`: the maximum battery charge power at time t.
- `Pbatt_dmax`: the maximum battery discharge power at time t.
- `Pgen_max`: the diesel generator rated power.
- `Pelyz_min` : the minimum electrolyzer power at time t.
- `Pelyz_max` : the maximum electrolyzer power at time t.
- 'Pfc_min' : the minimum electrolyzer power at time t.
- 'Pfc_max' : the minimum electrolyzer power at time t.
"""
function dispatch_2(Pnl_req,
                    Pbatt_cmax, Pbatt_dmax,
                    Pgen_min , Pgen_max,
                    Pelyz_min, Pelyz_max, 
                    Pfc_min,  Pfc_max)
    Pgen = 0.0
    Pspill = 0.0
    Pshed = 0.0
    Pelyz = 0.0
    Pfc = 0.0
    # Pnl_req >= 0 - load excess - after evaluating the production (Pnl = Pload - VRE generation)
    if Pnl_req >= 0.0
        # fuel_cell starts on if possible . the battery can't supply the power net load request
        if Pnl_req > Pfc_min
            Pfc = min(Pfc_max, Pnl_req)
        end

        # Battery starts on if the fuel cell can't supply the power net load request
        Pbatt = min(Pbatt_dmax, Pnl_req - Pfc)
   
        if Pnl_req - Pbatt - Pfc > Pgen_min         
            # Diesel generator starts on if the battery & the fuel_cell can't supply the power net load request
            Pgen = min(Pgen_max, Pnl_req - Pbatt - Pfc)
        end

        Pshed = Pnl_req - Pbatt - Pfc - Pgen

    # Pnl_req < 0 - VRE excess
    elseif Pnl_req < 0.0

        if -1 * Pnl_req > Pelyz_min
            #Electrolyzer starts on  
            Pelyz=min(Pelyz_max, - Pnl_req)
            
        end

        # Battery starts on if its remain enough power   
        Pbatt=max(Pbatt_cmax,Pnl_req + Pelyz)

        Pspill = Pbatt - Pelyz - Pnl_req
    end

    Pnl = Pnl_req - Pshed + Pspill

    return Pnl, Pgen , Pbatt, Pspill, Pshed, Pelyz, Pfc
end