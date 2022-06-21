"""
    operation(mg::Microgrid)

Simulate the annual operation of the microgrid `mg` and return the
hourly operation variables `OperVarsTraj`.
"""
function operation(mg::Microgrid)
    # photovoltaic production over one year
    # photovoltaic_production = production(mg.photovoltaic)
    renewables_production = collect(production(nd) for nd in mg.nondispatchables)

    # wind turbine production over one year
    # windpower_production = production(mg.windpower)

    # power balance before battery+generator
    # renewable_production = photovoltaic_production + windpower_production
    total_renewables_production = sum(renewables_production)::Vector{Float64}
    power_net_load_requested = mg.power_load - total_renewables_production

    # variables initialization
    stepsnumber = length(mg.power_load)
    T = Float64 # TODO: use some typeof call instead. But typeof which variable...?
    power_net_load = zeros(T,stepsnumber)
    Pgen = zeros(T,stepsnumber)
    Ebatt = zeros(T,stepsnumber+1)
    Ebatt[1] = mg.battery.energy_initial
    Pbatt = zeros(T,stepsnumber)
    Pbatt_dmax = zeros(T,stepsnumber)
    Pbatt_cmax = zeros(T,stepsnumber)
    Pcurt = zeros(T,stepsnumber)
    Pshed = zeros(T,stepsnumber)

    for i=1:stepsnumber
        # battery limits
        Pb_emin = - (mg.battery.energy_max - Ebatt[i]) / ((1 - mg.battery.loss) * mg.project.timestep)
        Pb_emax = (Ebatt[i] - mg.battery.energy_min) / ((1 + mg.battery.loss) * mg.project.timestep)
        Pbatt_dmax[i] = min(Pb_emax,mg.battery.power_max)
        Pbatt_cmax[i] = max(Pb_emin,mg.battery.power_min)
        
        # dispatch
        power_net_load[i], Pgen[i], Pbatt[i], Pcurt[i], Pshed[i] = dispatch(power_net_load_requested[i], Pbatt_cmax[i], Pbatt_dmax[i], mg.dieselgenerator.power_rated)

        # next battery energy level
        if i < stepsnumber
            Ebatt[i+1] = Ebatt[i] - (Pbatt[i] + mg.battery.loss * abs(Pbatt[i])) * mg.project.timestep
            # if Ebatt[i+1] < 0
            #     Ebatt[i+1] = 0
            # end
        end
    end
    
    opervarstraj = OperVarsTraj(power_net_load, Pshed, Pgen, Ebatt, Pbatt, Pbatt_dmax, Pbatt_cmax, Pcurt)
    return opervarstraj
end

"""
    aggregation(mg::Microgrid, opervarstraj::OperVarsTraj)

Return the aggregated operation variables `OperVarsAggr` for the microgrid `mg` and
the hourly operation variables `OperVarsTraj`.
"""
function aggregation(mg::Microgrid, opervarstraj::OperVarsTraj)
    # variables initialization
    nHours = 0
    vFuel  = 0
    Eserv = 0
    delestTimeCurrent = 0
    delestTimeMax = 0
    annualThrpt = 0.0
    for i = 1 : length(opervarstraj.Pgen)
        # diesel generator
        if opervarstraj.Pgen[i] > 0 #DG on
            nHours = nHours + 1 # incremention of the number of hours            
            fuelConsump = mg.dieselgenerator.F0 * mg.dieselgenerator.power_rated + mg.dieselgenerator.F1 * opervarstraj.Pgen[i] #[L.h-1]  F = F0*Ygen + F1*Pgen
        else
            fuelConsump = 0 # DG off [L.h-1]
        end
        vFuel = vFuel + fuelConsump*mg.project.timestep # [L]

        # battery
        annualThrpt = annualThrpt + abs(opervarstraj.Pbatt[i])*mg.project.timestep #battery throughtput (charging & dicharging)
        
        # load - energy served - energy provided (served in one year) [kWh]
        Eserv = Eserv + (mg.power_load[i] - opervarstraj.power_shedding[i]) * mg.project.timestep

        # load - calculation of the Max shedding duration
        if opervarstraj.power_shedding[i] > 0 && i!=length(opervarstraj.Pgen)
            delestTimeCurrent = delestTimeCurrent + 1
        else
            if delestTimeCurrent > delestTimeMax
                delestTimeMax = delestTimeCurrent
            end
            delestTimeCurrent = 0 # Ré-initialisation de la durée de délestage
        end        
    end

    # load
    Pshed_max = maximum(opervarstraj.power_shedding)    # max shedding
    shed_time_max = delestTimeMax
    shed_total_kwh = sum(opervarstraj.power_shedding) * mg.project.timestep # [kWh]
    shed_rate =  (sum(opervarstraj.power_shedding) * mg.project.timestep) / (sum(mg.power_load)*mg.project.timestep) # [0,1] 
    # battery
    annualThrpt = annualThrpt/2    # (charging and discharging cycle)
    # renewables sources
    Pcurt_max = maximum(opervarstraj.power_curtailment)    # max curtailment
    renew_rate = (1 - sum(opervarstraj.Pgen)/Eserv) * 100 # d'insertion ENR insertion Rate
    # outputs
    opervarsaggr = OperVarsAggr(Eserv, Pshed_max, shed_time_max, shed_total_kwh, shed_rate, nHours, vFuel, annualThrpt, Pcurt_max, renew_rate)
    
    # TODO - parei aqui
    # load
    # operVars.delestSumHours = sum(operVars.loaddelestRealHours); #[h]
    # operVars.delestRateHours =  100*sum(operVars.loaddelestRealHours) / length(mgDef.load); # [%]
    # operVars.serviceRate = 100 - operVars.load.delestRate; # service rate in term of Energy (%)
    # operVars.serviceRateHours = 100 - operVars.load.delestRateHours; #  service rate in term of duration[%]
    
    return opervarsaggr
end

"""
    aggregation(mg::Microgrid, opervarstraj::OperVarsTraj, ε)

When a relaxation factor `ε` is given, return the aggregated operation variables
`OperVarsAggr` with a linear relaxation of the diesel generator annual operation
hours `OperVarsAggr.nHours`.
"""
function aggregation(mg::Microgrid, opervarstraj::OperVarsTraj, ε)
    # variables initialization
    nHours = 0.
    vFuel  = 0
    Eserv = 0
    delestTimeCurrent = 0
    delestTimeMax = 0
    annualThrpt = 0
    for i = 1 : length(opervarstraj.Pgen)
        x = opervarstraj.Pgen[i]/mg.dieselgenerator.power_rated
        # diesel generator
        if opervarstraj.Pgen[i] > 0 && x <= ε   #DG on
            nHours = nHours + x/ε # incremention of the number of hours            
            fuelConsump = mg.dieselgenerator.F0 * mg.dieselgenerator.power_rated + mg.dieselgenerator.F1 * opervarstraj.Pgen[i] #[L.h-1]  F = F0*Ygen + F1*Pgen
        elseif opervarstraj.Pgen[i] > 0 && x > ε   #DG on
            nHours = nHours + 1 # incremention of the number of hours            
            fuelConsump = mg.dieselgenerator.F0 * mg.dieselgenerator.power_rated + mg.dieselgenerator.F1 * opervarstraj.Pgen[i] #[L.h-1]  F = F0*Ygen + F1*Pgen
        else
            fuelConsump = 0 # DG off [L.h-1]
        end
        vFuel = vFuel + fuelConsump*mg.project.timestep # [L]

        # battery
        annualThrpt = annualThrpt + abs(opervarstraj.Pbatt[i])*mg.project.timestep #battery throughtput (charging & dicharging)
        
        # load - energy served - energy provided (served in one year) [kWh]
        Eserv = Eserv + (mg.power_load[i] - opervarstraj.power_shedding[i]) * mg.project.timestep

        # load - calculation of the Max shedding duration
        if opervarstraj.power_shedding[i] > 0 && i!=length(opervarstraj.Pgen)
            delestTimeCurrent = delestTimeCurrent + 1
        else
            if delestTimeCurrent > delestTimeMax
                delestTimeMax = delestTimeCurrent
            end
            delestTimeCurrent = 0 # Ré-initialisation de la durée de délestage
        end        
    end

    # load
    Pshed_max = maximum(opervarstraj.power_shedding)    # max shedding
    shed_time_max = delestTimeMax
    shed_total_kwh = sum(opervarstraj.power_shedding) * mg.project.timestep # [kWh]
    shed_rate =  (sum(opervarstraj.power_shedding) * mg.project.timestep) / (sum(mg.power_load)*mg.project.timestep) * 100 # (%)
    # battery
    annualThrpt = annualThrpt/2    # (charging and discharging cycle)
    # renewables sources
    Pcurt_max = maximum(opervarstraj.power_curtailment)    # max curtailment
    renew_rate = (1 - sum(opervarstraj.Pgen .* mg.project.timestep)/Eserv) * 100 # d'insertion ENR insertion Rate
    # outputs
    opervarsaggr = OperVarsAggr(Eserv, Pshed_max, shed_time_max, shed_total_kwh, shed_rate, nHours, vFuel, annualThrpt, Pcurt_max, renew_rate)
    
    # TODO - parei aqui
    # load
    # operVars.delestSumHours = sum(operVars.loaddelestRealHours); #[h]
    # operVars.delestRateHours =  100*sum(operVars.loaddelestRealHours) / length(mgDef.load); # [%]
    # operVars.serviceRate = 100 - operVars.load.delestRate; # service rate in term of Energy (%)
    # operVars.serviceRateHours = 100 - operVars.load.delestRateHours; #  service rate in term of duration[%]
    
    return opervarsaggr
end

# # With the number of DG operation hours relaxed - exponential
# function aggregation(mg::Microgrid, opervarstraj::OperVarsTraj, ε)
#     # variables initialization
#     nHours = 0.
#     vFuel  = 0
#     Eserv = 0
#     delestTimeCurrent = 0
#     delestTimeMax = 0
#     annualThrpt = 0
#     for i = 1 : length(opervarstraj.Pgen)
#         x = opervarstraj.Pgen[i]/mg.dieselgenerator.power_rated
#         # diesel generator
#         if opervarstraj.Pgen[i] > 0   #DG on
#             nHours = nHours + exp(-ε * x) # incremention of the number of hours            
#             fuelConsump = mg.dieselgenerator.F0 * mg.dieselgenerator.power_rated + mg.dieselgenerator.F1 * opervarstraj.Pgen[i] #[L.h-1]  F = F0*Ygen + F1*Pgen
#         else
#             fuelConsump = 0 # DG off [L.h-1]
#         end
#         vFuel = vFuel + fuelConsump*mg.project.timestep # [L]

#         # battery
#         annualThrpt = annualThrpt + abs(opervarstraj.Pbatt[i])*mg.project.timestep #battery throughtput (charging & dicharging)
        
#         # load - energy served - energy provided (served in one year) [kWh]
#         Eserv = Eserv + (mg.load[i] - opervarstraj.Pshed[i]) * mg.project.timestep

#         # load - calculation of the Max shedding duration
#         if opervarstraj.Pshed[i] > 0 && i!=length(opervarstraj.Pgen)
#             delestTimeCurrent = delestTimeCurrent + 1
#         else
#             if delestTimeCurrent > delestTimeMax
#                 delestTimeMax = delestTimeCurrent
#             end
#             delestTimeCurrent = 0 # Ré-initialisation de la durée de délestage
#         end        
#     end

#     # load
#     Pshed_max = maximum(opervarstraj.Pshed)    # max shedding
#     shed_time_max = delestTimeMax
#     shed_total_kwh = sum(opervarstraj.Pshed) * mg.project.timestep # [kWh]
#     shed_rate =  (sum(opervarstraj.Pshed) * mg.project.timestep) / (sum(mg.load)*mg.project.timestep) # [0,1] 
#     # battery
#     annualThrpt = annualThrpt/2    # (charging and discharging cycle)
#     # renewables sources
#     Pcurt_max = maximum(opervarstraj.Pcurt)    # max curtailment
#     renew_rate = (1 - sum(opervarstraj.Pgen)/Eserv) * 100 # d'insertion ENR insertion Rate
#     # outputs
#     opervarsaggr = OperVarsAggr(Eserv, Pshed_max, shed_time_max, shed_total_kwh, shed_rate, nHours, vFuel, annualThrpt, Pcurt_max, renew_rate)
    
#     # TODO - parei aqui
#     # load
#     # operVars.delestSumHours = sum(operVars.loaddelestRealHours); #[h]
#     # operVars.delestRateHours =  100*sum(operVars.loaddelestRealHours) / length(mgDef.load); # [%]
#     # operVars.serviceRate = 100 - operVars.load.delestRate; # service rate in term of Energy (%)
#     # operVars.serviceRateHours = 100 - operVars.load.delestRateHours; #  service rate in term of duration[%]
    
#     return opervarsaggr
# end