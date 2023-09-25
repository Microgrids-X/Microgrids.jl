module Microgrids

import Base.@kwdef # backport Julia 1.9 syntax to 1.6-1.8 versions

export simulate,
       NonDispatchableSource,
       Project, DispatchableGenerator, Battery, Photovoltaic, PVInverter, WindPower, Microgrid,
       capacity_from_wind,
       OperationTraj, OperationStats,
       operation, aggregation, dispatch, production,
       CostFactors, MicrogridCosts, component_costs, economics

include("components.jl")
include("dispatch.jl")
include("production.jl")
include("operation.jl")
include("economics.jl")

function simulate(mg::Microgrid)
    # Run the microgrid operation
    opervarstraj = operation(mg)

    # Aggregate the operation variables
    opervarsaggr = aggregation(mg, opervarstraj)

    # Eval the microgrid costs
    costs = economics(mg, opervarsaggr)

    return (opervarstraj = opervarstraj, opervarsaggr = opervarsaggr, costs = costs)
end

end