module Microgrids

export simulate,
       NonDispatchables,
       Project, DieselGenerator, Battery, Photovoltaic, WindPower, Microgrid,
       TotalCosts, OperVarsTraj, OperVarsAggr,
       operation, aggregation, dispatch, production,
       annual_costs, economics

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