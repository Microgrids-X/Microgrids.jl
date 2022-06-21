using Microgrids
using BenchmarkTools
using CSV, DataFrames


function mg_create()
    data = DataFrame(CSV.File("$(@__DIR__)/../examples/microgrid_with_PV_BT_DG/data/Ouessant_data_2016.csv"))
    # Simulation steps
    ntimestep = length(data.Load)

    # Components parameters
    # Project
    lifetime = 25
    discount_rate = 0.05
    timestep = 1
    # Load
    Pload = data."Load"[1:ntimestep]
    # Photovoltaic
    power_rated_PV = 4106.82251423571
    fPV = 1.
    IT = data."Ppv1k"[1:ntimestep] ./ 1000
    IS = 1.
    investiment_cost_PV = 1200.
    om_cost_PV = 20.  
    replacement_cost_PV = 1200.
    salvage_cost_PV = 1200.
    lifetime_PV = 25
    # Battery
    energy_initial = 0.
    energy_max = 6839.87944197573
    energy_min = 0
    power_min = -1.114*energy_max
    power_max = 1.002*energy_max
    loss = 0.05
    investiment_cost_BT = 350.
    om_cost_BT = 10.
    replacement_cost_BT = 350.
    salvage_cost_BT = 350.
    lifetime_BT = 15
    lifetime_thrpt = 3000
    # Diesel generator
    power_rated_DG = 1800.
    min_load_ratio = 0
    F0 = 0.0
    F1 = 0.240
    fuel_cost = 1.
    investiment_cost_DG = 400.
    om_cost_DG = 0.02
    replacement_cost_DG = 400.
    salvage_cost_DG = 400.
    lifetime_DG = 15000

    # Create microgrid components
    project = Project(lifetime, discount_rate, timestep)
    dieselgenerator = DieselGenerator(power_rated_DG, min_load_ratio, F0, F1, fuel_cost, investiment_cost_DG, om_cost_DG, replacement_cost_DG, salvage_cost_DG, lifetime_DG)
    photovoltaic = Photovoltaic(power_rated_PV, fPV, IT, IS, investiment_cost_PV, om_cost_PV, replacement_cost_PV, salvage_cost_PV, lifetime_PV)
    battery = Battery(energy_initial, energy_max, energy_min, power_min, power_max, loss, investiment_cost_BT, om_cost_BT, replacement_cost_BT, salvage_cost_BT, lifetime_BT, lifetime_thrpt)

    # Create microgrid
    microgrid = Microgrid(project, Pload, dieselgenerator, battery, [photovoltaic])
end

mg = mg_create();
print("timing of simulate(mg):")
@btime results = simulate($mg)

function simulate_time(mg)
    # Run the microgrid operation
    opervarstraj = @btime operation($mg)

    # Aggregate the operation variables
    opervarsaggr = @btime aggregation($mg, $opervarstraj)

    # Eval the microgrid costs
    costs = @btime economics($mg, $opervarsaggr)

    return (opervarstraj = opervarstraj, opervarsaggr = opervarsaggr, costs = costs)
end

print("\ndetailed timing of simulate(mg):\n")
simulate_time(mg);