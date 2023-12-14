"""Operation variables (time series) from a simulated Microgrid operation

(simulation duration is assumed to be 1 year)
Works only for a `Microgrid` with 1 electrolyzer, one fuel_cell, and 1 one diesel generator .
"""
struct OperationTraj{T<:Real}
    # load
    #= "Net load at each time instant after using the renewables power (kW)"
    Pnl_req =#
    "Net load at each time instant after dispatch (kW)"
    power_net_load::Vector{T}
    "Unmet load/Load shedding power at each time instant (kW)"
    power_shedding::Vector{T}
    # non dispatchable (renewables)
    "renewable power potential (before spillage) (kW)"
    Prenew_pot::Vector{T}
    # dispatchables
    "Diesel generator power at each time instant (kW)"
    Pgen::Vector{T}
    "fuel cell power at each time instant (kW)"
    Pfc :: Vector{T}
     # electrolyzer
    "Electrolyzer power at each time instant (kW)"
    Pelyz:: Vector{T}
    # battery
    "Battery energy at each time instant (kWh)"
    Ebatt::Vector{T}
    "Battery power at each time instant (kW)"
    Pbatt::Vector{T}
    # tanks
    "Level of hydrogen at each time instant (kg)"
    LoH :: Vector{T}
    "Level of fuel at each time instant (L)"
    LoF :: Vector{T}
    # renewables sources
    "Renewables curtailment power at each time instant (kW)"
    power_curtailment::Vector{T}
end

"""Aggregated statistics over the simulated Microgrid operation

(simulation duration is assumed to be 1 year)
Works only for a `Microgrid` with 1 electrolyzer, one fuel_cell, and 1 one diesel generator .
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

    # Dispatchable diesel generator statistics
    "energy supplied by the dispatchable generator (kWh/y)"
    gen_energy
    "cumulated operating hours of the dispatchable generator (h/y)"
    gen_hours
    "fuel consumption (L/y)"
    gen_fuel

    # Hydrogen chain statistics
    "energy supplied by the fuel cell (kWh/y)"
    fc_energy
    "cumulated operating hours of the fuel_cell (h/y)"
    fc_hours
    "total amount of consumed hydrogen"
    h2_consumed
    "total amount of produced hydrogen"
    h2_produced
    "cumulated operating hours of the electrolyzer (h/y)"
    elyz_hours
    "energy consumed by the the electrolyzer to produce hydrogen (kWh/y)"
    elyz_consumed_energy
    "energy lossed in the Power-Gaz-Power chain (kWh/y)"
    h2_chain_loss

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

const oper_stats_units = Dict(
    :served_energy => "kWh",
    :shed_energy => "kWh",
    :shed_max => "kW",
    :shed_hours => "h",
    :shed_duration_max => "h",
    :shed_rate => "in [0,1]",
    :gen_energy => "kWh",
    :gen_hours => "h",
    :gen_fuel => "L",# TODO: use instead the generator's fuel unit string
    :fc_energy => "kWh",
    :fc_hours => "h",
    :h2_consumed => "kg",
    :h2_produced => "kg",
    :elyz_hours => "h",
    :elyz_consumed_energy => "kWh",
    :h2_chain_loss => "kWh",
    :storage_cycles => "",
    :storage_char_energy => "kWh",
    :storage_dis_energy => "kWh",
    :storage_loss_energy => "kWh",
    :spilled_energy => "kWh",
    :spilled_max => "kW",
    :spilled_rate => "in [0,1]",
    :renew_potential_energy => "kWh",
    :renew_energy => "kWh",
    :renew_rate => "in [0,1]"
)

"""
    Base.show(io::IO, ::MIME"text/plain", stats::OperationStats)

Write a multi-line text representation of `OperationStats` statistics
"""
function Base.show(io::IO, ::MIME"text/plain", stats::OperationStats)
    function format_stat(st::OperationStats, name::Symbol; sigdigits=6)
        value = getproperty(st, name)
        return round(value; sigdigits=sigdigits)
    end

    println(io, "OperationStats with fields:")
    for name in propertynames(stats)
        unit = get(oper_stats_units, name, "")
        println(io, "- ", name, ": ", format_stat(stats, name; sigdigits=5), " ", unit)
    end
end

"""
    increment(prod_unit_power::Float64, prod_unit::ProductionUnit,  ε::Real=0.0)

Allow to compute operating hours and fuel consumption of 'ProductionUnit'.
"""
function increment(prod_unit_power::Float64, prod_unit::ProductionUnit,  ε::Real=0.0)
    Pprod_unit_max = prod_unit.power_rated
    Pprod_unit_eps = Pprod_unit_max*1e-6
    prod_unit_intercept = prod_unit.consumption_intercept
    prod_unit_slope = prod_unit.consumption_slope
    time_inc = 0
    cons_rate = 0
        if prod_unit_power > Pprod_unit_eps # prod_unit ON
            Pprod_unit_norm = prod_unit_power / (ε * Pprod_unit_max) # can be Inf e.g. for ε=0.0
           
            if Pprod_unit_norm <= 1.0
                time_inc = Pprod_unit_norm # relaxation of discontinuous generator statistics for small Pgen
                cons_rate = Pprod_unit_norm * prod_unit_intercept * Pprod_unit_max +
                prod_unit_slope * prod_unit_power # (L/h)
            else
                time_inc=1
                cons_rate = prod_unit_intercept * Pprod_unit_max +
                prod_unit_slope * prod_unit_power # (L/h)
            end
          
        end
    return time_inc , cons_rate


end

"""
    operation(mg::Microgrid)

Simulate the annual operation of the microgrid `mg` and return the
hourly operation variables `OperationTraj`.
"""
function operation(mg::Microgrid)
    # Type of all variables: Float64 or ForwardDiff.Dual{...}
    Topt = typeof(mg).parameters[1]

    # Renewable power generation
    # (remark on naming convention: all non-dispatchable sources are assumed renewable!)
    renew_productions = collect(production(nd) for nd in mg.nondispatchables)
    renew_potential = sum(renew_productions)::Vector{Topt}

    # Desired net load
    Pnl_request = mg.load - renew_potential

    # Fixed parameters and short aliases
    K = length(mg.load)
    dt = mg.project.timestep

    Pgen_max = collect(gen.power_rated for gen in mg.dispatchables.generator)
    Pgen_min = collect(mg.dispatchables.generator[i].minimum_load_ratio * Pgen_max[i] for i in eachindex(Pgen_max) )
    fuel_slope = collect(gen.consumption_slope for gen in mg.dispatchables.generator)
    fuel_intercept = collect(gen.consumption_intercept for gen in mg.dispatchables.generator)

    LoF_max = mg.tanks.fuelTank.max_filling_ratio * mg.tanks.fuelTank.capacity
    LoF_min = mg.tanks.fuelTank.min_filling_ratio * mg.tanks.fuelTank.capacity

    Esto_max = mg.storage.energy_rated
    Esto_min = mg.storage.SoC_min * Esto_max
    Psto_pmax =  mg.storage.discharge_rate * Esto_max
    Psto_pmin = -mg.storage.charge_rate * Esto_max # <0 in line with the generator convention for Psto
    sto_loss = mg.storage.loss_factor

    n_elyz = length(mg.electrolyzer)
    Pelyz_max = collect(el.power_rated for el in mg.electrolyzer)
    Pelyz_min = collect(mg.electrolyzer[i].minimum_load_ratio * Pelyz_max[i] for i in eachindex(Pelyz_max) )
    prod_rate = collect(el.consumption_slope for el in mg.electrolyzer)

    Pfc_max = collect(fc.power_rated for fc in mg.dispatchables.fuel_cell)
    Pfc_min = collect(mg.dispatchables.fuel_cell[i].minimum_load_ratio * Pfc_max[i] for i in eachindex(Pfc_max) )
    cons_rate = collect(fc.consumption_slope for fc in mg.dispatchables.fuel_cell)

    LoH_max = mg.tanks.h2Tank.max_filling_ratio * mg.tanks.h2Tank.capacity
    LoH_min = mg.tanks.h2Tank.min_filling_ratio * mg.tanks.h2Tank.capacity

    

    # Initialization of loop variables
    Pnl = zeros(Topt,K)
    Pgen = zeros(Topt,K,length(mg.dispatchables.generator))
    

    Esto = zeros(Topt,K+1)
    Psto = zeros(Topt,K)



    Pspill = zeros(Topt,K)
    Pshed = zeros(Topt,K)

    Pelyz = zeros(Topt,K,length(mg.electrolyzer))
    Pfc = zeros(Topt,K,length(mg.dispatchables.fuel_cell))
    

    LoH = zeros(Topt,K+1)
    LoF = zeros(Topt,K+1)

    # Initial battery storage state
    Esto_ini = mg.storage.SoC_ini * mg.storage.energy_rated
    Esto[1] = Esto_ini

     # Initial h2 storage state
    LoH_ini = mg.tanks.h2Tank.ini_filling_ratio * mg.tanks.h2Tank.capacity
    LoH[1] = LoH_ini

    # Initial fuel tank state

    LoF_ini = mg.tanks.fuelTank.ini_filling_ratio * mg.tanks.fuelTank.capacity
    LoF[1] = LoF_ini
    
    Pelyz_emax = zeros(Topt, n_elyz) 

    for k=1:K
        # Storage energy and power limits 
        Psto_emin = - (Esto_max - Esto[k]) / ((1 - sto_loss) * dt)
        Psto_emax = (Esto[k] - Esto_min) / ((1 + sto_loss) * dt)
        Psto_dmax = min(Psto_emax, Psto_pmax)
        Psto_cmax = max(Psto_emin, Psto_pmin)


        #Pelyz_emax = collect( (LoH_max - LoH[k] ) * prod_rate[i] * dt for i in eachindex(prod_rate)) 
        Pelyz_emax .=  prod_rate .* ((LoH_max - LoH[k]) * dt)
        Pelyz_cmax =  collect( min(Pelyz_emax[i],Pelyz_max[i]) for  i in eachindex(Pelyz_emax))


        Pfc_emax = collect( (LoH[k] - LoH_min) / cons_rate[i] * dt for i in eachindex(cons_rate)) 
        Pfc_pmax =  collect( min(Pfc_emax[i],Pfc_max[i]) for  i in eachindex(Pfc_emax))

        Pgen_emax = collect( (((LoF[K] - LoF_min) - fuel_intercept[i] * Pgen_max[i]* dt )/ fuel_slope[i] * dt ) for i in eachindex(cons_rate)) 
        Pgen_dmax =  collect( min(Pgen_emax[i],Pgen_max[i]) for  i in eachindex(Pgen_emax))


        # dispatch
        Pnl[k], Pgen[k,1], Psto[k], Pspill[k], Pshed[k], Pelyz[k,1], Pfc[k,1] = dispatch(Pnl_request[k], Psto_cmax, Psto_dmax, Pgen_min[1], Pgen_dmax[1],  Pelyz_min[1],Pelyz_cmax[1], Pfc_min[1],Pfc_pmax[1])

        # Battery storage dynamics
        Esto[k+1] = Esto[k] - (Psto[k] + sto_loss * abs(Psto[k])) * dt

        # Hydrogen storage dynamics
        LoH[k+1] = LoH[k] + (
                sum(Pelyz[k,i]/prod_rate[i] for i in eachindex(prod_rate)) -
                sum(Pfc[k,i]*cons_rate[i] for i in eachindex(cons_rate))
            ) * dt

        # fuel storage dynamics
        LoF[k+1] = LoF[k] - sum(
                Pgen[k,i] / fuel_slope[i] + fuel_intercept[i]* Pgen_max[i]
                for i in eachindex(fuel_slope)
            ) * dt

    end

    oper_traj = OperationTraj(Pnl, Pshed, renew_potential, Pgen[:,1], Pfc[:,1], Pelyz[:,1], Esto, Psto, LoH, LoF, Pspill)
    
    return oper_traj
end

"""
aggregation(mg::Microgrid, oper_traj::OperationTraj, ε::Real=1.0)

Aggregate operation time series `oper_traj` into yearly statistics
for the the microgrid `mg` (returned as an `OperationStats` object).

Discontinuous statistics can optionally be relaxed (smoothed)
using the relaxation parameter `ε`:
- 0.0 means no relaxation (default value)
- 1.0 yields the strongest relaxation

Using relaxation (`ε` > 0) is recommended when using gradient-based optimization
and then a “small enough” value between 0.05 and 0.30 is suggested.
"""
function aggregation(mg::Microgrid, oper_traj::OperationTraj, ε::Real=0.0)
    ### Retrieve parameters
    dt = mg.project.timestep
    K = length(mg.load)

    ### Compute simple yearly statistics (sum, max...):
    # Load statistics
    load_energy = sum(mg.load) * dt # kWh/y
    shed_energy = sum(oper_traj.power_shedding) * dt # kWh/y
    served_energy = load_energy - shed_energy # kWh/y
    shed_max = maximum(oper_traj.power_shedding) # kW
    shed_rate = shed_energy / load_energy

    # Dispatchable generator statistics
    gen_energy = sum(oper_traj.Pgen) * dt # kWh/y

    # Energy storage (e.g. battery) statistics
    pos(x) = x >= 0.0 ? x : 0.0 # positive part
    neg(x) = x <= 0.0 ? -x : 0.0 # negative part
    storage_char_energy = sum(neg, oper_traj.Pbatt) * dt # kWh/y
    storage_dis_energy  = sum(pos, oper_traj.Pbatt) * dt # kWh/y
    storage_cycles = (storage_char_energy + storage_dis_energy) /
                     (2*mg.storage.energy_rated) # cycles/y
    Efin = oper_traj.Ebatt[end]
    Eini = oper_traj.Ebatt[1]
    storage_loss_energy = storage_char_energy - storage_dis_energy -
                          (Efin - Eini) # kWh/y

    #Power-H2-Power statistics  
    elyz_consumed_energy = sum(oper_traj.Pelyz) * dt
    elyz_produced_h2 = elyz_consumed_energy / mg.electrolyzer[1].consumption_slope
    fc_produced_energy = sum(oper_traj.Pfc) * dt
    fc_consumed_h2 = fc_produced_energy * mg.dispatchables.fuel_cell[1].consumption_slope
    
    LoHfin = oper_traj.LoH[end]
    LoHini = oper_traj.LoH[1]
    chain_loss =  elyz_consumed_energy - fc_produced_energy -(LoHfin-LoHini)*mg.dispatchables.fuel_cell[1].consumption_slope

    # Non-dispatchable (typ. renewables) sources statistics
    spilled_energy = sum(oper_traj.power_curtailment) * dt # kWh/y
    spilled_max = maximum(oper_traj.power_curtailment) # kW
    renew_potential_energy = sum(oper_traj.Prenew_pot) * dt # kWh/y
    spilled_rate = spilled_energy / renew_potential_energy
    renew_energy = renew_potential_energy - spilled_energy
    renew_rate = 1 - gen_energy/served_energy

    ### Iterative computation of more complex yearly statistics

    # Initialization of integrators:
    shed_hours = 0.0 # h/y
    shed_duration_max = 0.0 # h
    gen_hours = 0.0 # h/y
    gen_fuel = 0.0 # L/y
    elyz_hours = 0.0
    fc_hours = 0.0
   

    # auxilliary counter:
    shed_duration = 0.0; # duration of current load shedding event (h)
    Pload_avg = load_energy/(K*dt) # kW

    for K = 1:K
        # Dispatchable generator : operating hours and fuel consumption
        time_inc, cons_rate = increment(oper_traj.Pgen[K],mg.dispatchables.generator[1],ε) 
         gen_hours += time_inc * dt
         gen_fuel += cons_rate * dt


        time_inc, cons_rate = increment(oper_traj.Pelyz[K],mg.electrolyzer[1],ε)
        elyz_hours += time_inc * dt
        

        time_inc, cons_rate = increment(oper_traj.Pfc[K],mg.dispatchables.fuel_cell[1],ε)
        fc_hours += time_inc * dt
       
        # Load shedding: shedding duration and maximum shedding duration
        Pshed = oper_traj.power_shedding[K]
        if Pshed > 0.0
            Pshed_norm = Pshed / (ε * Pload_avg) # can be Inf e.g. for ε=0.0
            if Pshed_norm <= 1.0
                # relaxation of discontinuous load statistics for small Pshed
                shed_hours += Pshed_norm * dt
                shed_duration += Pshed_norm * dt
            else
                shed_hours += dt
                shed_duration += dt
            end
            # keep track of maximum shedding duration:
            shed_duration_max = max(shed_duration, shed_duration_max)
        else
            # reset shedding duration counter:
            shed_duration = 0.0
        end
    end

    # Outputs
    oper_stats = OperationStats(
        # Load statistics
        served_energy, shed_energy, shed_max, shed_hours, shed_duration_max, shed_rate,
        # Dispatchable generator statistics
        gen_energy, gen_hours, gen_fuel, fc_produced_energy, fc_hours, fc_consumed_h2,
        #
        elyz_produced_h2, elyz_hours, elyz_consumed_energy, chain_loss,

        # Energy storage (e.g. battery) statistics
        storage_cycles, storage_char_energy, storage_dis_energy, storage_loss_energy,
        # Non-dispatchable (typ. renewables) sources statistics
        spilled_energy, spilled_max, spilled_rate,
        renew_potential_energy, renew_energy, renew_rate
    )
    return oper_stats
end