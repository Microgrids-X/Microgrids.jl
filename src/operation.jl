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
    # non dispatchable (renewables)
    "renewable power potential (before spillage) (kW)"
    Prenew_pot::Vector{T}
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
    "fuel consumption (L/y)"
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

    oper_traj = OperationTraj(power_net_load, Pshed, Prenew_pot, Pgen, Ebatt, Pbatt, Pbatt_dmax, Pbatt_cmax, Pcurt)
    return oper_traj
end

"""
    aggregation(mg::Microgrid, oper_traj::OperationTraj, ε::Real=1.0)

Aggregates operation time series `oper_traj` into yearly statistics
for the the microgrid `mg` (returned as an `OperationStats` object).

Discontinuous statistics can optionally be relaxed (smoothed)
using the relaxation parameter `ε`:
- 0.0 means no relaxation (default value)
- 1.0 yields the strongest relaxation

when using relaxation, a value between 0.05 and 0.30 is suggested
"""
function aggregation(mg::Microgrid, oper_traj::OperationTraj, ε::Real=1.0)
    ### Retrieve parameters
    dt = mg.project.timestep
    nsteps = length(mg.load)

    ### Compute simple yearly statistics (sum, max...):
    # Load statistics
    load_energy = sum(mg.load) * dt # kWh/y
    shed_energy = sum(oper_traj.power_shedding) * dt # kWh/y
    served_energy = load_energy - shed_energy # kWh/y
    shed_max = max(oper_traj.power_shedding) # kW
    shed_rate = shed_energy / load_energy

    # Dispatchable generator statistics
    gen_energy = sum(oper_traj.Pgen) * dt # kWh/y

    # Energy storage (e.g. battery) statistics
    pos(x) = x >= 0.0 ? x : 0.0 # positive part
    storage_char_energy = sum(pos.(oper_traj.Pbatt)) * dt # kWh/y
    storage_dis_energy = sum(pos.(-oper_traj.Pbatt)) * dt # kWh/y
    storage_cycles = (storage_char_energy + storage_dis_energy) /
                     (2*mg.storage.energy_rated) # cycles/y
    Efin = oper_traj.Ebatt[end]
    Eini = oper_traj.Ebatt[1]
    storage_loss_energy = storage_char_energy - storage_dis_energy -
                          (Efin - Eini) # kWh/y

    # Non-dispatchable (typ. renewables) sources statistics
    spilled_energy = sum(oper_traj.power_curtailment) * dt # kWh/y
    spilled_max = max(oper_traj.power_curtailment) # kW
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

    # auxilliary counter:
    shed_duration = 0.0; # duration of current load shedding event (h)
    Pload_avg = load_energy/(nsteps*dt) # kW

    for i = 1:nsteps
        # Dispatchable generator : operating hours and fuel consumption
        Pgen = oper_traj.Pgen[i]
        Pgen_max = mg.generator.power_rated
        if Pgen > 0.0 # generator ON
            Pgen_norm = Pgen / (ε * Pgen_max) # can be Inf e.g. for ε=0.0
            if Pgen_norm <= 1.0
                # relaxation of discontinuous generator statistics for small Pgen
                gen_hours += Pgen_norm * dt
                fuel_rate = Pgen_norm * mg.generator.fuel_intercept * Pgen_max +
                            mg.generator.fuel_slope * Pgen # (L/h)
            else
                gen_hours += dt
                fuel_rate = mg.generator.fuel_intercept * Pgen_max +
                            mg.generator.fuel_slope * Pgen # (L/h)
            end
            gen_fuel += fuel_rate * dt
        end

        # Load shedding: shedding duration and maximum shedding duration
        Pshed = oper_traj.power_shedding[i]
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
        gen_energy, gen_hours, gen_fuel,
        # Energy storage (e.g. battery) statistics
        storage_cycles, storage_char_energy, storage_dis_energy, storage_loss_energy,
        # Non-dispatchable (typ. renewables) sources statistics
        spilled_energy, spilled_max, spilled_rate,
        renew_potential_energy, renew_energy, renew_rate
    )
    return oper_stats
end