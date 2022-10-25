using Microgrids
using BenchmarkTools
using CSV, DataFrames
using ForwardDiff

const data = DataFrame(CSV.File("$(@__DIR__)/../examples/microgrid_with_PV_BT_DG/data/Ouessant_data_2016.csv"))
# Simulation steps
const ntimestep = length(data.Load)

# Split Load and PV irradiance
const Pload = data."Load"[1:ntimestep]
const IT = data."Ppv1k"[1:ntimestep] ./ 1000

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
    power_rated_DG = x[1]
    energy_max = x[2]
    power_rated_PV = x[3]

    # Components parameters

    # Project
    lifetime = 25 # yr
    discount_rate = 0.05
    timestep = 1. # h

    # Diesel generator
    # power_rated_DG = 1800. # set in x
    min_load_ratio = 0.
    F0 = 0.0
    F1 = 0.240
    fuel_cost = 1.
    investiment_cost_DG = 400.
    om_cost_DG = 0.02
    replacement_cost_DG = 400.
    salvage_cost_DG = 400.
    lifetime_DG = 15000.

    # Battery energy storage
    energy_initial = 0.
    # energy_max = 6839.87944197573 # set in x
    energy_min = 0.
    power_min = -1.114*energy_max
    power_max = 1.002*energy_max
    loss = 0.05
    investiment_cost_BT = 350.
    om_cost_BT = 10.
    replacement_cost_BT = 350.
    salvage_cost_BT = 350.
    lifetime_BT = 15.
    lifetime_thrpt = 3000.

    # Photovoltaic generation
    # power_rated_PV = 4106.82251423571 # set in x
    fPV = 1.
    IS = 1.
    investiment_cost_PV = 1200.
    om_cost_PV = 20.
    replacement_cost_PV = 1200.
    salvage_cost_PV = 1200.
    lifetime_PV = 25.

    # Create microgrid components
    project = Project(lifetime, discount_rate, timestep)
    dieselgenerator = DieselGenerator(power_rated_DG, min_load_ratio, F0, F1, fuel_cost, investiment_cost_DG, om_cost_DG, replacement_cost_DG, salvage_cost_DG, lifetime_DG)
    photovoltaic = Photovoltaic(power_rated_PV, fPV, IT, IS, investiment_cost_PV, om_cost_PV, replacement_cost_PV, salvage_cost_PV, lifetime_PV)
    battery = Battery(energy_initial, energy_max, energy_min, power_min, power_max, loss, investiment_cost_BT, om_cost_BT, replacement_cost_BT, salvage_cost_BT, lifetime_BT, lifetime_thrpt)

    # Create microgrid
    microgrid = Microgrid(project, Pload, dieselgenerator, battery, [photovoltaic])
end

"""Microgrid with default sizing"""
mg_create() = mg_create([1800., 6839.9, 4106.8])

mg = mg_create();
print("timing of simulate(mg):")
@btime results = simulate($mg)

function simulate_time(mg)
    # Run the microgrid operation
    print("- operation:")
    opervarstraj = @btime operation($mg)

    # Aggregate the operation variables
    print("- aggregation:")
    opervarsaggr = @btime aggregation($mg, $opervarstraj)

    # Eval the microgrid costs
    print("- economics:")
    costs = @btime economics($mg, $opervarsaggr)

    return (opervarstraj = opervarstraj, opervarsaggr = opervarsaggr, costs = costs)
end

println("\ndetailed timing of simulate(mg):")
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