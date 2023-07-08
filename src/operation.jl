"""Operation variables (time series) from a simulated Microgrid operation

(simulation duration is assumed to be 1 year)
"""
struct OperationTraj{T<:Real}
    # load
    #= "Net load at each time instant after using the renewables power (kW)"
    Pnl_req =#
    "Net load at each time instant after dispatch (kW)"
    power_net_load::Vector{T}
    "Unmet load/Load shedding power at each time instant (kW)"
    power_shedding::Vector{T}
    # diesel generator
    "Diesel generator power at each time instant (kW)"
    Pgen::Vector{T}
    # battery
    "Battery energy at each time instant (kWh)"
    Ebatt::Vector{T}
    "Battery power at each time instant (kW)"
    Pbatt::Vector{T}
    "Maximum battery discharge power at time t (kW)"
    Pbatt_dmax::Vector{T}
    "Maximum battery charge power at time t (kW)"
    Pbatt_cmax::Vector{T}
    # renewables sources
    "Renewables curtailment power at each time instant (kW)"
    power_curtailment::Vector{T}
end

"""Aggregated statistics over the simulated Microgrid operation

(simulation duration is assumed to be 1 year)
"""
struct OperationStats
    # Load statistics
    "energy actually served to the load (kWh/y)"
    served_energy
    "shed energy, that is not served to the load (kWh/y)"
    shed_energy
    "maximum load shedding power (kW)"
    shed_max
    "cumulated duration of load shedding (h/y)"
    shed_hours
    "maximum consecutive duration of load shedding (h)"
    shed_duration_max
    "ratio of shed energy to the desired load (∈ [0,1])"
    shed_rate

    # Dispatchable generator statistics
    "energy supplied by the dispatchable generator (kWh/y)"
    gen_energy
    "cumulated operating hours of the dispatchable generator (h/y)"
    gen_hours
    "Diesel consumption in one year (L)"
    gen_fuel

    # Energy storage (e.g. battery) statistics
    "cycling of the energy storage (cycles/y)"
    storage_cycles
    "energy charged into the energy storage (kWh/y)"
    storage_char_energy
    "energy discharged out of the energy storage (kWh/y)"
    storage_dis_energy
    "energy lossed in the energy storage (kWh/y)"
    storage_loss_energy

    # Non-dispatchable (typ. renewables) sources statistics
    "spilled energy (typ. from excess of renewables) (kWh/y)"
    spilled_energy
    "maximum spilled power (typ. from excess of renewables) (kW)"
    spilled_max
    "ratio of spilled energy to the energy potentially supplied by renewables (∈ [0,1])"
    spilled_rate
    "energy potentially supplied by renewables in the absence spillage (kWh/y)"
    renew_potential_energy
    "energy actually supplied by renewables when substracting spillage (kWh/y)"
    renew_energy
    "ratio of energy actually supplied by renewables (net of storage loss) to the energy served to the load (∈ [0,1])"
    renew_rate
end


"""
    operation(mg::Microgrid)

Simulate the annual operation of the microgrid `mg` and return the
hourly operation variables `OperVarsTraj`.
"""
function operation(mg::Microgrid)
    # Type of all variables: Float64 or ForwardDiff.Dual{...}
    Topt = typeof(mg).parameters[1]

    # photovoltaic production over one year
    # photovoltaic_production = production(mg.photovoltaic)
    renewables_production = collect(production(nd) for nd in mg.nondispatchables)

    # wind turbine production over one year
    # windpower_production = production(mg.windpower)

    # power balance before battery+generator
    # renewable_production = photovoltaic_production + windpower_production

    total_renewables_production = sum(renewables_production)::Vector{Topt}
    power_net_load_requested = mg.power_load - total_renewables_production

    # variables initialization
    stepsnumber = length(mg.power_load)
    power_net_load = zeros(Topt,stepsnumber)
    Pgen = zeros(Topt,stepsnumber)
    Ebatt = zeros(Topt,stepsnumber+1)
    Ebatt[1] = mg.battery.energy_initial
    Pbatt = zeros(Topt,stepsnumber)
    Pbatt_dmax = zeros(Topt,stepsnumber)
    Pbatt_cmax = zeros(Topt,stepsnumber)
    Pcurt = zeros(Topt,stepsnumber)
    Pshed = zeros(Topt,stepsnumber)

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

    oper_traj = OperationTraj(power_net_load, Pshed, Pgen, Ebatt, Pbatt, Pbatt_dmax, Pbatt_cmax, Pcurt)
    return oper_traj
end

"""
    aggregation(mg::Microgrid, oper_traj::OperationTraj)

Return the aggregated operation statistics `OperationStats` for the microgrid `mg` and
the hourly operation variables `OperationTraj`.
"""
function aggregation(mg::Microgrid, oper_traj::OperationTraj)
    # variables initialization
    nHours = 0
    vFuel  = 0
    Eserv = 0
    delestTimeCurrent = 0
    delestTimeMax = 0
    annualThrpt = 0.0
    for i = 1 : length(oper_traj.Pgen)
        # diesel generator
        if oper_traj.Pgen[i] > 0 #DG on
            nHours = nHours + 1 # incremention of the number of hours
            fuelConsump = mg.dieselgenerator.F0 * mg.dieselgenerator.power_rated + mg.dieselgenerator.F1 * oper_traj.Pgen[i] #[L.h-1]  F = F0*Ygen + F1*Pgen
        else
            fuelConsump = 0 # DG off [L.h-1]
        end
        vFuel = vFuel + fuelConsump*mg.project.timestep # [L]

        # battery
        annualThrpt = annualThrpt + abs(oper_traj.Pbatt[i])*mg.project.timestep #battery throughtput (charging & dicharging)

        # load - energy served - energy provided (served in one year) [kWh]
        Eserv = Eserv + (mg.power_load[i] - oper_traj.power_shedding[i]) * mg.project.timestep

        # load - calculation of the Max shedding duration
        if oper_traj.power_shedding[i] > 0 && i!=length(oper_traj.Pgen)
            delestTimeCurrent = delestTimeCurrent + 1
        else
            if delestTimeCurrent > delestTimeMax
                delestTimeMax = delestTimeCurrent
            end
            delestTimeCurrent = 0 # Ré-initialisation de la durée de délestage
        end
    end

    # load
    Pshed_max = maximum(oper_traj.power_shedding)    # max shedding
    shed_time_max = delestTimeMax
    shed_total_kwh = sum(oper_traj.power_shedding) * mg.project.timestep # [kWh]
    shed_rate =  (sum(oper_traj.power_shedding) * mg.project.timestep) / (sum(mg.power_load)*mg.project.timestep) # [0,1]
    # battery
    annualThrpt = annualThrpt/2    # (charging and discharging cycle)
    # renewables sources
    Pcurt_max = maximum(oper_traj.power_curtailment)    # max curtailment
    renew_rate = (1 - sum(oper_traj.Pgen)/Eserv) * 100 # d'insertion ENR insertion Rate
    # outputs
    oper_stats = OperationStats(Eserv, Pshed_max, shed_time_max, shed_total_kwh, shed_rate, nHours, vFuel, annualThrpt, Pcurt_max, renew_rate)

    # TODO - parei aqui
    # load
    # operVars.delestSumHours = sum(operVars.loaddelestRealHours); #[h]
    # operVars.delestRateHours =  100*sum(operVars.loaddelestRealHours) / length(mgDef.load); # [%]
    # operVars.serviceRate = 100 - operVars.load.delestRate; # service rate in term of Energy (%)
    # operVars.serviceRateHours = 100 - operVars.load.delestRateHours; #  service rate in term of duration[%]

    return oper_stats
end

"""
    aggregation(mg::Microgrid, oper_traj::OperationTraj, ε)

When a relaxation factor `ε` is given, return the aggregated operation statistics
`OperationStats` with a linear relaxation of the diesel generator annual operation
hours `OperVarsAggr.nHours`.
"""
function aggregation(mg::Microgrid, oper_traj::OperationTraj, ε)
    # variables initialization
    nHours = 0.
    vFuel  = 0
    Eserv = 0
    delestTimeCurrent = 0
    delestTimeMax = 0
    annualThrpt = 0
    for i = 1 : length(oper_traj.Pgen)
        x = oper_traj.Pgen[i]/mg.dieselgenerator.power_rated
        # diesel generator
        if oper_traj.Pgen[i] > 0 && x <= ε   #DG on
            nHours = nHours + x/ε # incremention of the number of hours
            fuelConsump = mg.dieselgenerator.F0 * mg.dieselgenerator.power_rated + mg.dieselgenerator.F1 * oper_traj.Pgen[i] #[L.h-1]  F = F0*Ygen + F1*Pgen
        elseif oper_traj.Pgen[i] > 0 && x > ε   #DG on
            nHours = nHours + 1 # incremention of the number of hours
            fuelConsump = mg.dieselgenerator.F0 * mg.dieselgenerator.power_rated + mg.dieselgenerator.F1 * oper_traj.Pgen[i] #[L.h-1]  F = F0*Ygen + F1*Pgen
        else
            fuelConsump = 0 # DG off [L.h-1]
        end
        vFuel = vFuel + fuelConsump*mg.project.timestep # [L]

        # battery
        annualThrpt = annualThrpt + abs(oper_traj.Pbatt[i])*mg.project.timestep #battery throughtput (charging & dicharging)

        # load - energy served - energy provided (served in one year) [kWh]
        Eserv = Eserv + (mg.power_load[i] - oper_traj.power_shedding[i]) * mg.project.timestep

        # load - calculation of the Max shedding duration
        if oper_traj.power_shedding[i] > 0 && i!=length(oper_traj.Pgen)
            delestTimeCurrent = delestTimeCurrent + 1
        else
            if delestTimeCurrent > delestTimeMax
                delestTimeMax = delestTimeCurrent
            end
            delestTimeCurrent = 0 # Ré-initialisation de la durée de délestage
        end
    end

    # load
    Pshed_max = maximum(oper_traj.power_shedding)    # max shedding
    shed_time_max = delestTimeMax
    shed_total_kwh = sum(oper_traj.power_shedding) * mg.project.timestep # [kWh]
    shed_rate =  (sum(oper_traj.power_shedding) * mg.project.timestep) / (sum(mg.power_load)*mg.project.timestep) * 100 # (%)
    # battery
    annualThrpt = annualThrpt/2    # (charging and discharging cycle)
    # renewables sources
    Pcurt_max = maximum(oper_traj.power_curtailment)    # max curtailment
    renew_rate = (1 - sum(oper_traj.Pgen .* mg.project.timestep)/Eserv) * 100 # d'insertion ENR insertion Rate
    # outputs
    oper_stats = OperationStats(Eserv, Pshed_max, shed_time_max, shed_total_kwh, shed_rate, nHours, vFuel, annualThrpt, Pcurt_max, renew_rate)

    # TODO - parei aqui
    # load
    # operVars.delestSumHours = sum(operVars.loaddelestRealHours); #[h]
    # operVars.delestRateHours =  100*sum(operVars.loaddelestRealHours) / length(mgDef.load); # [%]
    # operVars.serviceRate = 100 - operVars.load.delestRate; # service rate in term of Energy (%)
    # operVars.serviceRateHours = 100 - operVars.load.delestRateHours; #  service rate in term of duration[%]

    return oper_stats
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
