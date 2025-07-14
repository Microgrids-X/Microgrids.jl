module Microgrids

import Base.@kwdef # backport Julia 1.9 syntax to 1.6-1.8 versions

export simulate,
       SalvageType, LinearSalvage, ConsistentSalvage,
       Project, Microgrid, Component,
       DispatchableGenerator, Battery,
       NonDispatchableSource, Photovoltaic, PVInverter, WindPower,
       capacity_from_wind,
       Smoothing, NoSmoothing, OperationTraj, OperationStats,
       operation, aggregation, dispatch, production,
       CostFactors, MicrogridCosts, component_costs, economics

include("components.jl")
include("dispatch.jl")
include("production.jl")
include("operation.jl")
include("economics.jl")

"""
    simulate(mg::Microgrid, smoothing::Smoothing=NoSmoothing)

Simulate the technical and economic performance of a Microgrid `mg`.

Discontinuous computations can optionally be smoothed (relaxed) using the
`Smoothing` parameters `transition` and `gain` (see their doc).
Smoothing is recommended when using gradient-based optimization.

Returns:
- Operational trajectories from `operation` (should be optional in future version)
- Operational statistics from `operation`
- Microgrid project costs from `economics`
"""
function simulate(mg::Microgrid, smoothing::Smoothing=NoSmoothing)
    # Run the microgrid operation
    oper_traj = operation(mg)

    # Aggregate the operation variables
    oper_stats = aggregation(mg, oper_traj, smoothing)

    # Eval the microgrid costs
    mg_costs = economics(mg, oper_stats)

    return (traj=oper_traj, stats=oper_stats, costs=mg_costs)
end

end