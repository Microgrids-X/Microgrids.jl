"""Recorder for trajectories of operational variables"""
struct TrajRecorder{T<:Real}
    values::Dict{Symbol,Vector{T}}
end

"""Initialize a TrajRecorder for trajectories of operational variables
of given names, with arrays prescribed length.

Example: recorder = TrajRecorder(var1=10, var2=11)
"""
function TrajRecorder(T; kwargs...)
    recorder = TrajRecorder{T}(Dict{Symbol,Vector{T}}())

    for (name, traj_length) in kwargs
        recorder.values[name] = zeros(T, traj_length)
    end

    return recorder
end

"""Record trajectory values at index `k`.

`k` should be in [1, `length`], for the given `length` provided at initialization.

Example: rec(recorder, 5, var1=0.8)
"""
function rec(recorder::TrajRecorder, k::Int; kwargs...)
    for (name, value) in kwargs
        recorder.values[name][k] = value
    end
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
hourly operation variables `OperVarsTraj`. TODO: update doctring.
"""
function operation(mg::Microgrid, ε::Real=0.0, record=false::Bool)
    # Type of all variables: Float64 or ForwardDiff.Dual{...}
    Topt = typeof(mg).parameters[1]

    # Renewable power generation
    # (remark on naming convention: all non-dispatchable sources are assumed renewable!)
    renew_productions = collect(production(nd) for nd in mg.nondispatchables)
    renew_potential = sum(renew_productions)::Vector{Topt}

    # TODO future version: no need to create renew_potential array
    # (including allocating each production(nd) array)

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

    if record
        recorder = TrajRecorder(Topt; Prep=K, Pgen=K, Psto=K, Esto=K+1, Pspill=K, Pshed=K)
    end

    for k=1:K
        ### Decide energy dispatch
        # Desired net load
        Pnl_request = mg.load[k] - renew_potential[k]

        # Storage energy and power limits
        Psto_emin = - (Esto_max - Esto) / ((1 - sto_loss) * dt)
        Psto_emax = (Esto - Esto_min) / ((1 + sto_loss) * dt)
        Psto_dmax = min(Psto_emax, Psto_pmax)
        Psto_cmax = max(Psto_emin, Psto_pmin)

        # Dispatch
        Pnl, Pgen, Psto, Pspill, Pshed = dispatch(Pnl_request, Psto_cmax, Psto_dmax, Pgen_max)

        # Record trajectories (optional)
        if record
            rec(recorder, k,
                Prep=renew_potential[k],
                Pgen=Pgen, Psto=Psto, Esto=Esto,
                Pspill=Pspill, Pshed=Pshed)
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

    if record
        rec(recorder, K+1, Esto=Esto) # Esto at last instant
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