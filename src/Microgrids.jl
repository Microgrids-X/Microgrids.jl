module Microgrids
using CSV, DataFrames
using DocStringExtensions



@template (TYPES) =
"""
$(TYPEDEF)
$(DOCSTRING)

### Parameters: 
$(FIELDS)
"""

include("components.jl")
include("dispatch.jl")
include("production.jl")
include("operation.jl")
include("economics.jl")
include("Microgrid_Wind-Solar_data.jl")

import Base.@kwdef # backport Julia 1.9 syntax to 1.6-1.8 versions

export simulate,
       NonDispatchableSource, ProductionUnit, Tank, TankCompound, 
       Project, DispatchableCompound,Battery, Photovoltaic, PVInverter, WindPower, Microgrid,
       capacity_from_wind,
       OperationTraj, OperationStats, increment,
       operation, aggregation, dispatch_1,dispatch_2, production, nsteps,
       CostFactors, MicrogridCosts, component_costs, economics, new_microgrid,
       energy_rated_sto, power_rated_pv,power_rated_wind, power_rated_elyz,power_rated_fc ,capacity_rated_hytank,Pload,X,capex_def



"""
    simulate(mg::Microgrid, ε::Real=0.0)

Simulate the technical and economic performance of a Microgrid `mg`.

Discontinuous computations can optionally be relaxed (smoothed)
using the relaxation parameter `ε`:
- 0.0 means no relaxation (default value)
- 1.0 yields the strongest relaxation

Using relaxation (`ε` > 0) is recommended when using gradient-based optimization
and then a “small enough” value between 0.05 and 0.30 is suggested.

Returns:
- Operational trajectories from `sim_operation` (should be optional in future version)
- Operational statistics from `sim_operation`
- Microgrid project costs from `sim_economics`
"""
function simulate(mg::Microgrid, dispatch::Function,ε::Real=0.0)
    # Run the microgrid operation
    oper_traj = operation(mg, dispatch)
    # Aggregate the operation variables
    oper_stats = aggregation(mg, oper_traj, ε)

    # Eval the microgrid costs
    mg_costs = economics(mg, oper_stats)

    return oper_traj, oper_stats, mg_costs
end



  """
    new_microgrid(x::Vector{Float64}= X,capex::Vector{Float64}=capex_def,initial_fill_rate::Vector{Float64})

TBW
"""
function new_microgrid(x::Vector{Float64}= X,capex::Vector{Float64}=capex_def,initial_fill_rate::Vector{Float64}=ini_filling_state)
    project = Project(lifetime, discount_rate, timestep, "€")
  gen = ProductionUnit(power_rated_gen,
      fuel_intercept, fuel_slope, fuel_price,
      capex[1], om_price_gen, lifetime_gen,
      load_ratio_min_gen, replacement_price_ratio,
      salvage_price_ratio, input_unit_gen,output_unit_gen)
  ftank = Tank(capacity_rated_ftank,capex[2], om_price_ftank,lifetime_ftank,loss_factor_ftank,initial_fill_rate[1],
          fuel_min_ratio, fuel_max_ratio,fuel_price, replacement_price_ratio, salvage_price_ratio)
  fuel_cell = ProductionUnit(x[5]*1000,cons_intercept_fc, cons_rate_fc,cons_price_fc,capex[3], om_price_fc,lifetime_fc,
              load_min_ratio_fc,replacement_price_ratio, salvage_price_ratio,input_unit_fc,output_unit_fc)
  hytank = Tank(x[6]*1000,capex[4], om_price_hytank,lifetime_hytank,loss_factor_hytank,initial_fill_rate[2],
              LoH_min_ratio, LoH_max_ratio,hy_price,replacement_price_ratio, salvage_price_ratio)
  dispatchables = DispatchableCompound{Float64}([gen], [fuel_cell])
     tanks = TankCompound{Float64}(ftank,hytank)

  elyz = ProductionUnit(x[4]*1000,cons_intercept_elyz,cons_slope_elyz,cons_price_elyz, capex[5], om_price_elyz, lifetime_elyz,
  load_min_ratio_elyz,replacement_price_ratio, salvage_price_ratio,input_unit_elyz,output_unit_elyz)

  batt = Battery(x[1]*1000,
      capex[6], om_price_sto, lifetime_sto, lifetime_cycles,
      charge_rate, discharge_rate, loss_factor_sto, SoC_min, SoC_max,initial_fill_rate[3],
      replacement_price_ratio, salvage_price_ratio)
  pv = Photovoltaic(x[2]*1000, irradiance,
      capex[7], om_price_pv,
      lifetime_pv, derating_factor_pv,
      replacement_price_ratio, salvage_price_ratio)
  windgen = WindPower(x[3]*1000, cf_wind,
      capex[8], om_price_wind,
      lifetime_wind,
      replacement_price_ratio, salvage_price_ratio)
  
 mg = Microgrid(project, Pload,dispatchables,
  [elyz,],tanks,batt, [
      pv,
      windgen
      ])
  return mg
  end
end