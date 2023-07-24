using Microgrids
using BenchmarkTools
using CSV, DataFrames
using ForwardDiff

const data = DataFrame(CSV.File("$(@__DIR__)/../examples/microgrid_with_PV_BT_DG/data/Ouessant_data_2016.csv"))
# Simulation steps
const nsteps = length(data.Load)

# Split Load and PV irradiance
const Pload = data."Load"[1:nsteps]
const irradiance = data."Ppv1k"[1:nsteps] ./ 1000

"""
    mg_create(x)

Create microgrid of size `x` (size=3).
Size vector `x` should contain:
- power_rated_DG
- energy_max (battery capacity)
- power_rated_PV
"""
function mg_create(x)
    # Split sizing vector:
    power_rated_gen = x[1]
    energy_rated = x[2]
    power_rated_pv = x[3]

    # Components parameters

    # Project
    lifetime = 25 # yr
    discount_rate = 0.05
    timestep = 1. # h

    # Parameters common to all Components
    replacement_price_ratio = 1.0
    salvage_price_ratio = 1.0

    # Dispatchable generator
    # power_rated_gen = 1800. # set in x
    fuel_intercept = 0.0
    fuel_slope = 0.240
    fuel_price = 1.
    investment_price_gen = 400.
    om_price_gen = 0.02
    lifetime_gen = 15000.
    load_ratio_min = 0.0
    fuel_unit = "L"

    # Battery energy storage
    # energy_rated = 6839.9 # set in x
    investment_price_sto = 350.
    om_price_sto = 10.
    lifetime_sto = 15.
    lifetime_cycles = 3000.
    charge_rate = 1.0
    discharge_rate = 1.0
    loss_factor_sto = 0.05
    SoC_min = 0.
    SoC_ini = 0.

    # Photovoltaic generation
    # power_rated_pv = 4106.8 # set in x
    investment_price_pv = 1200.
    om_price_pv = 20.
    lifetime_pv = 25.
    derating_factor_pv = 1.

    # Components:
    generator = DispatchableGenerator(power_rated_gen,
        fuel_intercept, fuel_slope, fuel_price,
        investment_price_gen, om_price_gen, lifetime_gen,
        load_ratio_min,
        replacement_price_ratio, salvage_price_ratio, fuel_unit)
    battery = Battery(energy_rated,
        investment_price_sto, om_price_sto, lifetime_sto, lifetime_cycles,
        charge_rate, discharge_rate, loss_factor_sto, SoC_min, SoC_ini,
        replacement_price_ratio, salvage_price_ratio)
    photovoltaic = Photovoltaic(power_rated_pv, irradiance,
        investment_price_pv, om_price_pv,
        lifetime_pv, derating_factor_pv,
        replacement_price_ratio, salvage_price_ratio)

    # Microgrid project, with discount rate:
    proj = Project(lifetime, discount_rate, timestep, "€")
    mg = Microgrid(proj, Pload, generator, battery, [photovoltaic])

    return mg
end

"""Microgrid with default sizing"""
mg_create() = mg_create([1800., 6839.9, 4106.8])

mg = mg_create();
print("timing of simulate(mg):")
@btime results = simulate($mg)

function simulate_time(mg)
    # Run the microgrid operation (including aggregation)
    print("- operation:")
    oper_stats = @btime operation($mg)

    # Eval the microgrid costs
    print("- economics:")
    costs = @btime economics($mg, $oper_stats)

    return (oper_stats = oper_stats, costs = costs)
end

#println("\ndetailed timing of simulate(mg):")
simulate_time(mg);

## Now performance of differentiation ##

"""simulate microgrid of size `x` (size=3) and returns its Net Present Cost"""
function sim_npc(x)
    mg = mg_create(x)

    results = simulate(mg)
    return results.costs.npc
end

print("timing of sim_npc(x):")
@btime sim_npc([1800., 6839.9, 4106.8]) # 183 µs, 167 alloc

print("timing of gradient(sim_npc, x):")
@btime ForwardDiff.gradient(sim_npc, [1800., 6839.9, 4106.8]) # 350 μs, i.e. 1.9×Tsim (181 allocations: 2.75 MiB)