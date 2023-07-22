module Microgrids

import Base.@kwdef # backport Julia 1.9 syntax to 1.6-1.8 versions

export simulate,
       NonDispatchables,
       Project, DispatchableGenerator, Battery, Photovoltaic, PVInverter, WindPower, Microgrid,
       OperationTraj, OperationStats,
       operation, aggregation, dispatch, production,
       CostFactors, MicrogridCosts, component_costs, economics

include("components.jl")
include("dispatch.jl")
include("production.jl")
include("operation.jl")
include("economics.jl")

function simulate(mg::Microgrid)
    # Run the microgrid operation and returs aggregated stats
    oper_stats = operation(mg)

    # Eval the microgrid costs
    costs = economics(mg, oper_stats)

    return (oper_stats = oper_stats, costs = costs)
end

end