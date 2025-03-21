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
    "haber_bosch power at each time instant (kW)"
    Phb:: Vector{T}
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
    Pdump::Vector{T}
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
    "cumulated starts on of the dispatchable generator "
    gen_starts
    "fuel consumption (L/y)"
    gen_fuel
    "running time at nominal power"
    gen_rt_at_nom
    "minimum running power"
    gen_min_rp
    "average running power"
    gen_avg_rp
    "average running time by per start"
    gen_avg_rt
    "maximum running time per start"
    gen_max_rt
    "minimum running time per start"
    gen_min_rt
    "number of runs that last 2h or less"
    gen_ru2
    "number of runs that last more than 2h"
    gen_ro2 

    # Hydrogen chain statistics
    "energy supplied by the fuel cell (kWh/y)"
    fc_energy
    "cumulated operating hours of the fuel_cell (h/y)"
    fc_hours
    "number of starts"
    fc_starts
    "total amount of consumed hydrogen"
    h2_consumed
    "running time at nominal power"
    fc_rt_at_nom
    "minimum running power"
    fc_min_rp
    "average running power"
    fc_avg_rp
    "average running time per start"
    fc_avg_rt
    "maximum running time per start"
    fc_max_rt
    "minimum running time per start"
    fc_min_rt
    "number of runs that last 2h or less"
    fc_ru2
    "number of runs that last more than 2h"
    fc_ro2
    "total amount of produced hydrogen"
    h2_produced
    "cumulated operating hours of the electrolyzer (h/y)"
    elyz_hours
    "number of starts"
    elyz_starts
    "energy consumed by the the electrolyzer to produce hydrogen (kWh/y)"
    elyz_consumed_energy
    "running time at nominal power"
    elyz_rt_at_nom
    "minimum running power"
    elyz_min_rp
    "average running power"
    elyz_avg_rp
    "average running time per start"
    elyz_avg_rt
    "maximum running time per start"
    elyz_max_rt
    "minimum running time per start"
    elyz_min_rt
    "number of runs that last 2h or less"
    elyz_ru2
    "number of runs that last more than 2h"
    elyz_ro2
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
    
    # haber-bosch statistics
    "energy supplied to the hb (kWh/y)"
    hb_cons_el
    "hydrogen supplied to the hb (kg/y)"
    hb_cons_h2
    "cumulated operating hours of the dispatchable generator (h/y)"
    hb_hours
    "cumulated starts on of the dispatchable generator "
    hb_starts
    "fuel consumption (L/y)"
    hb_prod
    "running time at nominal power"
    hb_rt_at_nom
    "minimum running power"
    hb_min_rp
    "average running power"
    hb_avg_rp
    "average running time by per start"
    hb_avg_rt
    "maximum running time per start"
    hb_max_rt
    "minimum running time per start"
    hb_min_rt
    "number of runs that last 2h or less"
    hb_ru2
    "number of runs that last more than 2h"
    hb_ro2 

    # Non-dispatchable (typ. renewables) sources statistics
    "spilled energy (typ. from excess of renewables) (kWh/y)"
    spilled_energy
    "maximum spilled power (typ. from excess of renewables) (kW)"
    spilled_max
    "ratio of spilled energy to the energy potentially supplied by renewables (∈ [0,1])"
    spilled_rate
    "energy potentially supplied by renewables in the absence spillage (kWh/y)"
    dumped_energy
    dump_max
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
    :gen_starts => "",
    :gen_fuel => "L",# TODO: use instead the generator's fuel unit string
    :gen_rt_at_nom =>"h",
    :gen_min_rp =>"KW",
    :gen_avg_rp => "kW",
    :gen_avg_rt => "h",
    :gen_max_rt =>"h",
    :gen_min_rt =>"h",
    :gen_ru2 =>"",
    :gen_ro2 =>"",
    :fc_energy => "kWh",
    :fc_hours => "h",
    :fc_starts => "",
    :h2_consumed => "kg",
    :fc_rt_at_nom =>"h",
    :fc_min_rp =>"KW",
    :fc_avg_rp => "kW",
    :fc_avg_rt => "h",
    :fc_max_rt =>"h",
    :fc_min_rt =>"h",
    :fc_ru2 =>"",
    :fc_ro2 =>"",
    :h2_produced => "kg",
    :elyz_hours => "h",
    :elyz_starts => "",
    :elyz_consumed_energy => "kWh",
    :elyz_rt_at_nom =>"h",
    :elyz_min_rp =>"KW",
    :elyz_avg_rp => "kW",
    :elyz_avg_rt => "h",
    :elyz_max_rt =>"h",
    :elyz_min_rt =>"h",
    :elyz_ru2 =>"",
    :elyz_ro2 =>"",
    :h2_chain_loss => "kWh",
    :storage_cycles => "",
    :storage_char_energy => "kWh",
    :storage_dis_energy => "kWh",
    :storage_loss_energy => "kWh",
    :hb_cons_el => "kWh",
    :hb_cons_h2 => "kWh",
    :hb_hours => "h",
    :hb_starts => "",
    :hb_prod => "kg",
    :hb_rt_at_nom => "h",
    :hb_min_rp => "h",
    :hb_avg_rp => "kWh",
    :hb_avg_rt => "h",
    :hb_max_rt => "h",
    :hb_min_rt => "h",
    :hb_ru2 => "",
    :hb_ro2  => "",
    :spilled_energy => "kWh",
    :spilled_max => "kW",
    :spilled_rate => "in [0,1]",
    :dump_energy => "kWh",
    :dump_max => "kW",
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
    usage(prod_unit_power::Float64,prod_unit::ProductionUnit,hours::Float64,cons::Float64,minP::Float64,dt::Float64,c1::Float64,max_rt::Float64,min_rt::Float64,r_u_2::Int64,r_o_2::Int64,cr::Float64,n_starts::Int64,ε::Real=0.0)

Allow to operate running statistics of production units
"""
function usage!(prod_unit_power::Float64,prod_unit::ProductionUnit,hours::Float64,cons::Float64,minP::Float64,dt::Float64,c1::Float64,max_rt::Float64,min_rt::Float64,r_u_2::Int64,r_o_2::Int64,cr::Float64,n_starts::Int64,ε::Real=0.0)
    Pprod_unit_max = prod_unit.power_rated
    Pprod_unit_eps = Pprod_unit_max*1e-5
    prod_unit_intercept = prod_unit.consumption_intercept
    prod_unit_slope = prod_unit.consumption_slope
    time_inc = 0.0
    cons_rate = 0
    if prod_unit_power > Pprod_unit_eps # prod_unit ON
        Pprod_unit_norm = prod_unit_power / (ε * Pprod_unit_max) # can be Inf e.g. for ε=0.0
        if Pprod_unit_norm <= 1.0
            time_inc = Pprod_unit_norm*dt # relaxation of discontinuous prod_unit statistics for small Prod_unit_power
            cons_rate = Pprod_unit_norm * prod_unit_intercept * Pprod_unit_max +
            prod_unit_slope * prod_unit_power # (L/h)
        else
            time_inc=dt
            cons_rate = prod_unit_intercept * Pprod_unit_max +
            prod_unit_slope * prod_unit_power # (L/h)
        end
        if cr == 0.0
            n_starts+=1
        end
        if prod_unit_power >= Pprod_unit_max - Pprod_unit_eps
            c1+=dt
        end
        if prod_unit_power < minP
            minP=prod_unit_power
        end
        cr+=dt
    else
        if cr>max_rt
            max_rt=cr
        end
        if cr<min_rt && cr>0.0
            min_rt=cr
        end
        if cr<2*dt && cr>0.0
            r_u_2+=1
        elseif cr>=2*dt && cr>0.0
            r_o_2+=1
        end
        cr=0
    end
    hours+=time_inc*dt
    cons+=cons_rate*dt
    return n_starts,minP,c1,max_rt,min_rt,r_u_2,r_o_2,cr,hours,cons
end
"""
    operation(mg::Microgrid, dispatch::Function)

Simulate the annual operation of the microgrid `mg` and return the
hourly operation variables `OperationTraj`.
"""
function operation(mg::Microgrid, dispatch::Function)
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

    n_gen = length(mg.dispatchables.generator)
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

    n_fc = length(mg.dispatchables.fuel_cell)
    Pfc_max = collect(fc.power_rated for fc in mg.dispatchables.fuel_cell)
    Pfc_min = collect(mg.dispatchables.fuel_cell[i].minimum_load_ratio * Pfc_max[i] for i in eachindex(Pfc_max) )
    cons_rate = collect(fc.consumption_slope for fc in mg.dispatchables.fuel_cell)

    LoH_max = mg.tanks.h2Tank.max_filling_ratio * mg.tanks.h2Tank.capacity
    LoH_min = mg.tanks.h2Tank.min_filling_ratio * mg.tanks.h2Tank.capacity

    # Initialization of loop variables
    Pnl = zeros(Topt,K)
    
    Pgen = zeros(Topt,K,n_elyz)
    Pgen_dmax = zeros(Topt, n_gen) 
    Pgen_emax = zeros(Topt, n_gen)
    
    Esto = zeros(Topt,K+1)
    Psto = zeros(Topt,K)

    Pspill = zeros(Topt,K)
    Pshed = zeros(Topt,K)
    Pdump = zeros(Topt,K)

    Pelyz = zeros(Topt,K, n_elyz)
    Pelyz_emax = zeros(Topt, n_elyz) 
    Pelyz_cmax = zeros(Topt, n_elyz)

    Pfc = zeros(Topt,K, n_fc)
    Pfc_emax = zeros(Topt, n_fc) 
    Pfc_pmax = zeros(Topt, n_fc)

    Phb = zeros(Topt,K)
    LoH = zeros(Topt,K+1)
    LoF = zeros(Topt,K+1)

    # Initial battery storage state
    Esto_ini = mg.storage.SoC_ini * mg.storage.energy_rated
    Esto[1] = Esto_ini

    # Initial H2 storage state
    LoH_ini = mg.tanks.h2Tank.ini_filling_ratio * mg.tanks.h2Tank.capacity
    LoH[1] = LoH_ini

    # Initial fuel tank state
    LoF_ini = mg.tanks.fuelTank.ini_filling_ratio * mg.tanks.fuelTank.capacity
    LoF[1] = LoF_ini

   
  
    for k=1:K
        # Storage energy and power limits 
        Psto_emin = - (Esto_max - Esto[k]) / ((1 - sto_loss) * dt)
        Psto_emax = (Esto[k] - Esto_min) / ((1 + sto_loss) * dt)
        Psto_dmax = min(Psto_emax, Psto_pmax)
        Psto_cmax = max(Psto_emin, Psto_pmin)
        
        # Electrolyzers power limits
        Pelyz_emax .=  prod_rate .* ((LoH_max - LoH[k]) * dt)
        Pelyz_cmax .=  min.(Pelyz_emax,Pelyz_max)

        # fuel_cells power limits
        Pfc_emax .=  max.(dt * (LoH[k] - LoH_min) ./ cons_rate , 0.0 ) 
        Pfc_pmax .=  min.(Pfc_emax,Pfc_max) 

        # Diesel generators power limits 
        Pgen_emax .= ((LoF[K] - LoF_min)/dt .- fuel_intercept.* Pgen_max )./ (fuel_slope)
        Pgen_dmax .=  min.(Pgen_emax,Pgen_max)

        # dispatch
        Pnl[k], Pgen[k,1], Psto[k], Pspill[k], Pshed[k], Pdump[k], Pelyz[k,1], Pfc[k,1] = dispatch(Pnl_request[k], Psto_cmax, Psto_dmax, Pgen_min[1], Pgen_dmax[1],  Pelyz_min[1],Pelyz_cmax[1], Pfc_min[1],Pfc_pmax[1])

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

    oper_traj = OperationTraj(Pnl, Pshed, renew_potential, Pgen[:,1], Pfc[:,1], Pelyz[:,1],Phb[:,1], Esto, Psto, LoH, LoF, Pspill,Pdump)
    
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
    co= 8760/K

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

    # Power-H2-Power statistics  
    elyz_consumed_energy = sum(oper_traj.Pelyz) * dt
    elyz_produced_h2 = elyz_consumed_energy / mg.electrolyzer[1].consumption_slope
    fc_produced_energy = sum(oper_traj.Pfc) * dt
    fc_consumed_h2 = fc_produced_energy * mg.dispatchables.fuel_cell[1].consumption_slope
    
    LoHfin = oper_traj.LoH[end]
    LoHini = oper_traj.LoH[1]
    chain_loss =  elyz_consumed_energy - fc_produced_energy -(LoHfin-LoHini)*mg.dispatchables.fuel_cell[1].consumption_slope
     # Haber-bosch statistics 
     hb_cons_el = sum(oper_traj.Phb) * dt
     hb_prod = hb_cons_el / mg.haber_bosch.consumption_slope
     hb_cons_h2=  hb_prod*0.176
     
    # Non-dispatchable (typ. renewables) sources statistics
    spilled_energy = sum(oper_traj.power_curtailment) * dt # kWh/y
    spilled_max = maximum(oper_traj.power_curtailment) # kW
    dumped_energy = sum(oper_traj.Pdump)*dt  # kWh/y
    dumped_max = maximum(oper_traj.Pdump)# kW
    renew_potential_energy = sum(oper_traj.Prenew_pot) * dt # kWh/y
    spilled_rate = spilled_energy / renew_potential_energy
    renew_energy = renew_potential_energy - spilled_energy
    renew_rate = 1 - gen_energy/served_energy

    ### Iterative computation of more complex yearly statistics

    # Initialization of integrators:
    shed_hours = 0.0 # h/y
    shed_duration_max = 0.0 # h
    gen_hours = 0.0 # h/y
    gen_starts = 0 # 
    gen_fuel = 0.0 # L/y
    gen_rt_at_nom = 0.0;
    gen_min_rp= mg.dispatchables.generator[1].power_rated;
    gen_avg_rp= 0.0;
    gen_avg_rt=0.0;
    gen_max_rt=0.0;
    gen_min_rt=K*dt
    gen_ru2=0;
    gen_ro2=0;
    elyz_hours = 0.0
    elyz_starts=0
    elyz_cons=0.0;
    elyz_rt_at_nom=0.0;
    elyz_min_rp=mg.electrolyzer[1].power_rated;
    elyz_avg_rp=0.0;
    elyz_avg_rt=0.0;
    elyz_max_rt=0.0;
    elyz_min_rt=K*dt
    elyz_ru2=0;
    elyz_ro2=0;
    fc_hours = 0.0
    fc_starts=0
    fc_cons=0.0;
    fc_rt_at_nom=0.0;
    fc_min_rp=mg.dispatchables.fuel_cell[1].power_rated;
    fc_avg_rp=0.0;
    fc_avg_rt=0.0;
    fc_max_rt=0.0;
    fc_min_rt=K*dt;
    fc_ru2=0;
    fc_ro2=0;
    hb_hours=0.;
    hb_starts=0;
    hb_rt_at_nom=0.;
    hb_min_rp=mg.haber_bosch.power_rated;
    hb_avg_rp=0.;
    hb_avg_rt=0.;
    hb_max_rt=0.;
    hb_min_rt=0.;
    hb_ru2=0;
    hb_ro2 =0;

    # auxilliary counter:
    shed_duration = 0.0; # duration of current load shedding event (h)
    Pload_avg = load_energy/(K*dt) # kW
    cr=zeros(Float64,4)
    for k = 1:K
        gen_starts,gen_min_rp,gen_rt_at_nom,gen_max_rt,gen_min_rt,gen_ru2,gen_ro2,cr[1],gen_hours,gen_fuel=usage!(oper_traj.Pgen[k],mg.dispatchables.generator[1],gen_hours,gen_fuel,gen_min_rp,dt,
                                                                                                                     gen_rt_at_nom,gen_max_rt,gen_min_rt,gen_ru2,gen_ro2,cr[1],gen_starts)

        fc_starts,fc_min_rp,fc_rt_at_nom,fc_max_rt,fc_min_rt,fc_ru2,fc_ro2,cr[2],fc_hours,fc_cons=usage!(oper_traj.Pfc[k],mg.dispatchables.fuel_cell[1],fc_hours,fc_cons,fc_min_rp,dt,
                                                                                                                    fc_rt_at_nom,fc_max_rt,fc_min_rt,fc_ru2,fc_ro2,cr[2],fc_starts)
        elyz_starts,elyz_min_rp,elyz_rt_at_nom,elyz_max_rt,elyz_min_rt,elyz_ru2,elyz_ro2,cr[3],elyz_hours,elyz_cons=usage!(oper_traj.Pelyz[k],mg.electrolyzer[1],elyz_hours,elyz_cons,elyz_min_rp,dt,
                                                                                                                    elyz_rt_at_nom,elyz_max_rt,elyz_min_rt,elyz_ru2,elyz_ro2,cr[3],elyz_starts)

        hb_starts,hb_min_rp,hb_rt_at_nom,hb_max_rt,hb_min_rt,hb_ru2,hb_ro2,cr[4],hb_hours,hb_cons_el=usage!(oper_traj.Phb[k],mg.haber_bosch,hb_hours,hb_cons_el,hb_min_rp,dt,
                                                                                                                    hb_rt_at_nom,hb_max_rt,hb_min_rt,hb_ru2,hb_ro2,cr[4],hb_starts)
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

    gen_avg_rp= sum(oper_traj.Pgen)/gen_hours;
    gen_avg_rt=gen_hours/gen_starts;
    fc_avg_rp= sum(oper_traj.Pfc)/fc_hours;
    fc_avg_rt=fc_hours/fc_starts;
    elyz_avg_rp= sum(oper_traj.Pelyz)/elyz_hours;
    elyz_avg_rt=elyz_hours/elyz_starts;
    hb_avg_rp= sum(oper_traj.Phb)/hb_hours;
    hb_avg_rt= hb_hours/hb_starts;

    # Outputs
oper_stats = OperationStats(
        # Load statistics
        served_energy*co, shed_energy*co, shed_max, shed_hours*co, shed_duration_max, shed_rate,
        # Dispatchable generator statistics
        gen_energy*co, gen_hours*co, gen_starts*co,gen_fuel*co,gen_rt_at_nom*co,gen_min_rp,gen_avg_rp,gen_avg_rt,gen_max_rt,gen_min_rt,gen_ru2*co,gen_ro2*co, 
        #H2 chain statistics
        fc_produced_energy*co, fc_hours*co,fc_starts*co, fc_consumed_h2*co,fc_rt_at_nom*co,fc_min_rp,fc_avg_rp,fc_avg_rt,fc_max_rt,fc_min_rt,fc_ru2*co,fc_ro2*co,
        elyz_produced_h2*co, elyz_hours*co,elyz_starts*co, elyz_consumed_energy*co,elyz_rt_at_nom*co,elyz_min_rp,elyz_avg_rp,elyz_avg_rt,elyz_max_rt,elyz_min_rt,elyz_ru2*co,elyz_ro2*co,
        chain_loss*co,

        # Energy storage (e.g. battery) statistics
        storage_cycles*co, storage_char_energy*co, storage_dis_energy*co, storage_loss_energy*co,
         # Energy storage (e.g. battery) statistics
        hb_cons_el, hb_cons_h2,hb_hours,hb_starts,hb_prod,hb_rt_at_nom,hb_min_rp,hb_avg_rp,hb_avg_rt,
        hb_max_rt,hb_min_rt,hb_ru2,hb_ro2,

        # Non-dispatchable (typ. renewables) sources statistics
        spilled_energy*co, spilled_max, spilled_rate,dumped_energy*co,dumped_max,
        renew_potential_energy*co, renew_energy*co, renew_rate
    )
    return oper_stats
end
