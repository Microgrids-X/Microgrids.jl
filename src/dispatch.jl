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
    Pdump = 0.0
    
    if Pnl_req >= 0.0
        ## Pnl_req >= 0 - load excess - after evaluating the production (Pnl = Pload - VRE generation)
        # battery discharging --> Pbatt positive

        if Pnl_req > Pbatt_dmax && Pfc_max >= Pfc_min
            if Pnl_req - Pbatt_dmax >= Pfc_min
                Pbatt =  Pbatt_dmax
                Pfc = min(Pfc_max,Pnl_req-Pbatt)
                Pshed = Pnl_req - Pbatt - Pfc - Pgen 
            else 
                Pfc=Pfc_min
                if Pnl_req>=Pfc
                    Pbatt= Pnl_req-Pfc
                    Pshed = 0.0
                else
                    
               Pbatt = max(Pnl_req-Pfc,Pbatt_cmax)
                Pdump = Pnl_req-Pfc-Pbatt
                end
            end
        else
            Pbatt=min(Pnl_req,Pbatt_dmax)
            Pshed = Pnl_req - Pbatt - Pfc - Pgen 
        end

    """    if Pnl_req - Pbatt - Pfc > Pgen_min
            # Diesel generator starts on if the battery & the fuel_cell can't supply the power net load request
            Pgen = min(Pgen_max,Pnl_req-Pbatt-Pfc)
        end
        """
 
    

    ## Pnl_req < 0 - VRE excess
    elseif Pnl_req < 0.0
        # battery charging --> Pbatt negative
        Pbatt = max(Pbatt_cmax,Pnl_req)
       
        if Pbatt - Pnl_req > Pelyz_min && Pelyz_max > Pelyz_min
            #Electrolyzer starts on if its remain enough power    
            Pelyz = min(Pelyz_max,Pbatt-Pnl_req)
        end

        Pspill = Pbatt - Pelyz - Pnl_req
    end

    Pnl = Pnl_req - Pshed + Pspill

    return Pnl, Pgen , Pbatt, Pspill, Pshed,Pdump, Pelyz, Pfc
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
    Pdump = 0.0

    # Pnl_req >= 0 - load excess - after evaluating the production (Pnl = Pload - VRE generation)
    if Pnl_req >= 0.0
        if Pfc_max >= Pfc_min
            if Pnl_req >= Pfc_min
                Pfc=min(Pfc_max,Pnl_req)
                Pbatt=min(Pnl_req-Pfc,Pbatt_dmax)
            else 
                if Pbatt_dmax < Pnl_req
                    Pfc=Pfc_min
                    Pbatt=max(Pnl_req-Pfc,Pbatt_cmax)
                    Pdump = Pnl_req-Pfc-Pbatt
                else
                    Pbatt = Pnl_req
                end

            end
        else
            Pbatt=min(Pnl_req,Pbatt_dmax)
        end
        
"""
        if Pnl_req - Pbatt - Pfc > Pgen_min
            # Diesel generator starts on if the battery & the fuel_cell can't supply the power net load request
            Pgen = min(Pgen_max,Pnl_req-Pbatt-Pfc)
        end
"""
        Pshed = Pnl_req - Pbatt - Pfc - Pgen - Pdump
    # Pnl_req < 0 - VRE excess
    elseif Pnl_req < 0.0

        if -1 * Pnl_req > Pelyz_min && Pelyz_max >= Pelyz_min
            #Electrolyzer starts on  
            Pelyz=min(Pelyz_max, -1* Pnl_req)
            
        end

        # Battery starts on if its remain enough power   
        Pbatt=max(Pbatt_cmax,Pnl_req + Pelyz)

        Pspill = Pbatt - Pelyz - Pnl_req
    end

    Pnl = Pnl_req - Pshed + Pspill

    return Pnl, Pgen , Pbatt, Pspill, Pshed, Pdump,Pelyz, Pfc
end