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
@kwdef mutable struct OperationStats{T<:Real}
    # Load statistics
    "energy actually served to the load (kWh/y)"
    served_energy::T
    "shed energy, that is not served to the load (kWh/y)"
    shed_energy::T
    "maximum load shedding power (kW)"
    shed_max::T
    "cumulated duration of load shedding (h/y)"
    shed_hours::T
    "maximum consecutive duration of load shedding (h)"
    shed_duration_max::T
    "ratio of shed energy to the desired load (∈ [0,1])"
    shed_rate::T

    # Dispatchable generator statistics
    "energy supplied by the dispatchable generator (kWh/y)"
    gen_energy::T
    "cumulated operating hours of the dispatchable generator (h/y)"
    gen_hours::T
    "fuel consumption (L/y)"
    gen_fuel::T

    # Energy storage (e.g. battery) statistics
    "cycling of the energy storage (cycles/y)"
    storage_cycles::T
    "energy charged into the energy storage (kWh/y)"
    storage_char_energy::T
    "energy discharged out of the energy storage (kWh/y)"
    storage_dis_energy::T
    "energy lossed in the energy storage (kWh/y)"
    storage_loss_energy::T

    # Non-dispatchable (typ. renewables) sources statistics
    "spilled energy (typ. from excess of renewables) (kWh/y)"
    spilled_energy::T
    "maximum spilled power (typ. from excess of renewables) (kW)"
    spilled_max::T
    "ratio of spilled energy to the energy potentially supplied by renewables (∈ [0,1])"
    spilled_rate::T
    "energy potentially supplied by renewables in the absence spillage (kWh/y)"
    renew_potential_energy::T
    "energy actually supplied by renewables when substracting spillage (kWh/y)"
    renew_energy::T
    "ratio of energy actually supplied by renewables (net of storage loss) to the energy served to the load (∈ [0,1])"
    renew_rate::T
end


"""
    operation(mg::Microgrid)

Simulate the annual operation of the microgrid `mg` and return the
hourly operation variables `OperVarsTraj`.
"""
function operation(mg::Microgrid, ε::Real=0.0, recorder=false)
    # Type of all variables: Float64 or ForwardDiff.Dual{...}
    Topt = typeof(mg).parameters[1]

    # Renewable power generation
    # (remark on naming convention: all non-dispatchable sources are assumed renewable!)
    renew_productions = collect(production(nd) for nd in mg.nondispatchables)
    renew_potential = sum(renew_productions)::Vector{Topt}

    # Desired net load
    Pnl_request = mg.load - renew_potential

    # TODO future version : no need to create renew_potential and Pnl_request arrays

    # Fixed parameters and short aliases
    K = length(mg.load)
    dt = mg.project.timestep
    Pgen_max = mg.generator.power_rated
    Esto_max = mg.storage.energy_rated
    Esto_min = mg.storage.SoC_min * Esto_max
    Psto_pmax =  mg.storage.discharge_rate * Esto_max
    Psto_pmin = -mg.storage.charge_rate * Esto_max # <0 in line with the generator convention for Psto
    sto_loss = mg.storage.loss_factor

    # Load statistics
    load_energy = sum(mg.load)*dt
    Pload_avg = load_energy/(K*dt) # kW

    # Initialization of loop variables
    # Initial storage state
    Esto_ini = mg.storage.SoC_ini * mg.storage.energy_rated
    Esto = Esto_ini
    # Operation statistics counters initialiazed at zero
    op_st = OperationStats{Topt}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0,0.0, 0.0)
    shed_duration = 0.0 # duration of current load shedding event (h)

    if recorder
        #init
    end

    for k=1:K
        ### Decide energy dispatch
        # Storage energy and power limits
        Psto_emin = - (Esto_max - Esto) / ((1 - sto_loss) * dt)
        Psto_emax = (Esto - Esto_min) / ((1 + sto_loss) * dt)
        Psto_dmax = min(Psto_emax, Psto_pmax)
        Psto_cmax = max(Psto_emin, Psto_pmin)

        # dispatch
        Pnl, Pgen, Psto, Pspill, Pshed = dispatch(Pnl_request[k], Psto_cmax, Psto_dmax, Pgen_max)

        # optionally record values in trajectories
        if recorder
            #record
        end

        # Storage dynamics
        Esto = Esto - (Psto + sto_loss*abs(Psto)) * dt

        ### Aggregate operation statistics
        # Load statistics
        if Pshed > 0.0 # load shedding
            op_st.shed_energy += Pshed*dt
            op_st.shed_max = max(op_st.shed_max, Pshed)

            if ε > 0.0 # relaxation of discontinuities ON
                Pshed_norm = Pshed / (ε * Pload_avg)
                if Pshed_norm <= 1.0
                    # relaxation of discontinuous load statistics for small Pshed
                    op_st.shed_hours += Pshed_norm * dt
                    shed_duration += Pshed_norm * dt
                else
                    op_st.shed_hours += dt
                    shed_duration += Pshed_norm * dt
                end
            else # relaxation OFF
                op_st.shed_hours += dt
                shed_duration += Pshed_norm * dt
            end

            op_st.shed_duration_max = max(op_st.shed_duration_max, shed_duration)

        else # no load shedding
            # reset duration counter of current load shedding event:
            shed_duration = 0.0
        end

        # Dispatchable generator statistics
        if Pgen > 0.0 # Generator ON
            op_st.gen_energy += Pgen*dt
            if ε > 0.0 # relaxation of discontinuities ON
                Pgen_norm = Pgen / (ε * Pgen_max)
                if Pgen_norm <= 1.0
                    # relaxation of discontinuous generator statistics for small Pgen
                    op_st.gen_hours += Pgen_norm * dt
                    fuel_rate = Pgen_norm * mg.generator.fuel_intercept * Pgen_max +
                                mg.generator.fuel_slope * Pgen # (L/h)
                else
                    op_st.gen_hours += dt
                    fuel_rate = mg.generator.fuel_intercept * Pgen_max +
                                  mg.generator.fuel_slope * Pgen # (L/h)
                end
            else # relaxation OFF
                op_st.gen_hours += dt
                fuel_rate = mg.generator.fuel_intercept * Pgen_max +
                              mg.generator.fuel_slope * Pgen # (L/h)
            end
            op_st.gen_fuel += fuel_rate*dt
        end

        # Energy storage (e.g. battery) statistics
        if Psto > 0.0 # discharge
            op_st.storage_dis_energy += Psto*dt
        else
            op_st.storage_char_energy -= Psto*dt
        end

        # Non-dispatchable (typ. renewables) sources statistics
        op_st.spilled_energy += Pspill*dt
        op_st.spilled_max = max(op_st.spilled_max, Pspill)

    end # for each instant k

    if recorder
        # record last instant of Esto (length K+1)
    end

    # Some more aggregated operation statistics
    op_st.served_energy = load_energy - op_st.shed_energy
    op_st.shed_rate = op_st.shed_energy / load_energy

    op_st.storage_loss_energy = op_st.storage_char_energy -
        op_st.storage_dis_energy - (Esto - Esto_ini)
    storage_throughput = op_st.storage_char_energy + op_st.storage_dis_energy
    op_st.storage_cycles = storage_throughput / (2*mg.storage.energy_rated)

    op_st.renew_potential_energy = sum(renew_potential)
    op_st.renew_energy = op_st.renew_potential_energy - op_st.spilled_energy
    op_st.renew_rate = 1 - op_st.gen_energy/op_st.served_energy
    op_st.spilled_rate = op_st.spilled_energy / op_st.renew_potential_energy

    return op_st
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

    # auxilliary counter:
    shed_duration = 0.0; # duration of current load shedding event (h)
    Pload_avg = load_energy/(K*dt) # kW

    for k = 1:K
        # Dispatchable generator : operating hours and fuel consumption
        Pgen = oper_traj.Pgen[k]
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
        Pshed = oper_traj.power_shedding[k]
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