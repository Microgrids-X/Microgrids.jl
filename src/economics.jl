# Economic modeling of a microgrid project

# TODO 2022: use the ComponentCosts struct as an output of annual_costs functions, rather than Vectors
"Cost factors of a Microgrid component, expressed as Net Present Values over the Microgrid project lifetime"
struct ComponentCosts
    "Total cost (initial + replacement + O&M + fuel + salvage)"
    total
    "Initial investment cost"
    investment
    "Replacement cost"
    replacement
    "Operation & Maintenance (O&M) cost"
    om
    "Fuel cost"
    fuel
    "Salvage cost (negative)"
    salvage
end

# TODO 2022: split the giant MicrogridCosts struct into a hierarchical struct of structs
"Cost components of a Microgrid project"
struct MicrogridCosts
    # general
    "Levelized cost of electricity (currency unit)"
    lcoe
    "Cost of electricity (currency unit)"
    coe # annualized
    "Net present cost (currency unit)"
    npc
    "Present investment cost (currency unit)"
    total_investment_cost
    "Present replacement cost (currency unit)"
    total_replacement_cost
    "Present operation and maintenance cost (currency unit)"
    total_om_cost
    "Present salvage cost (currency unit)"
    total_salvage_cost

    # components
    "Generator's total present cost (currency unit)"
    DG_total_cost
    "Generator's present investment cost (currency unit)"
    DG_investment_cost
    "Generator's present replacement cost (currency unit)"
    DG_replacement_cost
    "Generator's present operation and maintenance cost (currency unit)"
    DG_om_cost
    "Generator's present salvage cost (currency unit)"
    DG_salvage_cost
    "Generator's present fuel cost (currency unit)"
    DG_fuel_cost

    "Battery's total present cost (currency unit)"
    BT_total_cost
    "Battery's present investment cost (currency unit)"
    BT_investment_cost
    "Battery's present replacement cost (currency unit)"
    BT_replacement_cost
    "Battery's present operation and maintenance cost (currency unit)"
    BT_om_cost
    "Battery's present salvage cost (currency unit)"
    BT_salvage_cost

    "Photovoltaic's total present cost (currency unit)"
    PV_total_cost
    "Photovoltaic's present investment cost (currency unit)"
    PV_investment_cost
    "Photovoltaic's present replacement cost (currency unit)"
    PV_replacement_cost
    "Photovoltaic's present operation and maintenance cost (currency unit)"
    PV_om_cost
    "Photovoltaic's present salvage cost (currency unit)"
    PV_salvage_cost

    "Wind turbine's total present cost (currency unit)"
    WT_total_cost
    "Wind turbine's present investment cost (currency unit)"
    WT_investment_cost
    "Wind turbine's present replacement cost (currency unit)"
    WT_replacement_cost
    "Wind turbine's present operation and maintenance cost (currency unit)"
    WT_om_cost
    "Wind turbine's present salvage cost (currency unit)"
    WT_salvage_cost
end


function annual_costs(mg_project::Project, quantity, investment_price, replacement_price, salvage_price, om_price, fuel_consumption, fuel_price, lifetime)
    # discount factor for each year of the project
    discount_factors = [ 1/((1 + mg_project.discount_rate)^i) for i=1:mg_project.lifetime ]
    sum_discounts = sum(discount_factors)

    # number of replacements
    replacements_number = ceil(Integer, mg_project.lifetime/lifetime) - 1
    # years that the replacements happen
    replacement_years = [i*lifetime for i=1:replacements_number]
    # discount factors for the replacements years
    replacement_factors = [1/(1 + mg_project.discount_rate)^i for i in replacement_years]

    # component remaining life at the project end
    remaining_life = lifetime*(1+replacements_number) - mg_project.lifetime
    # proportional unitary salvage cost given remaining life
    salvage_price_effective = salvage_price * remaining_life / lifetime

    # present investment cost
    investment_cost = investment_price * quantity
    # present operation and maintenance cost
    om_cost = om_price * quantity * sum_discounts
    # present replacement cost
    if replacements_number == 0
        replacement_cost = 0.0
    else
        replacement_cost = replacement_price * quantity * sum(replacement_factors)
    end
    # Salvage cost (<0)
    salvage_cost = -salvage_price_effective * quantity * discount_factors[mg_project.lifetime]

    if fuel_consumption > 0.0
        fuel_cost = fuel_price * fuel_consumption * sum_discounts
    else
        fuel_cost = 0.0
    end

    total_cost = investment_cost + replacement_cost + om_cost + fuel_cost + salvage_cost

    return ComponentCosts(total_cost, investment_cost, replacement_cost, om_cost, fuel_cost, salvage_cost)
end

"""costs for NonDispatchables (PV, wind...) components"""
function annual_costs(nd::NonDispatchables, mg_project::Project)
    c = annual_costs(
        mg_project,
        nd.power_rated,
        nd.investment_cost,
        nd.replacement_cost,
        nd.salvage_cost,
        nd.om_cost,
        0.0, 0.0,
        nd.lifetime)
    return [c.total, c.investment, c.om, c.replacement, -c.salvage]
end

function annual_costs(pvi::PVInverter, mg_project::Project)
    c_ac = annual_costs(
        mg_project,
        pvi.power_rated,
        pvi.investment_cost_ac,
        pvi.replacement_cost_ac,
        pvi.salvage_cost_ac,
        pvi.om_cost_ac,
        0.0, 0.0,
        pvi.lifetime_ac)
    c_dc = annual_costs(
        mg_project,
        pvi.power_rated*pvi.ILR, #DC Power Rated
        pvi.investment_cost_dc,
        pvi.replacement_cost_dc,
        pvi.salvage_cost_dc,
        pvi.om_cost_dc,
        0.0, 0.0,
        pvi.lifetime_dc)
    return [c_ac.total + c_dc.total,
            c_ac.investment + c_dc.investment,
            c_ac.om + c_dc.om,
            c_ac.replacement+c_dc.replacement,
            -(c_ac.salvage+c_dc.salvage)]
end

function annual_costs(dg::DieselGenerator, mg_project::Project, opervarsaggr::OperVarsAggr)

    # discount factor for each year of the project
    discount_factors = [ 1/((1 + mg_project.discount_rate)^i) for i=1:mg_project.lifetime ]

    # total diesel generator operation hours over the project lifetime
    total_DG_operation_hours = mg_project.lifetime * opervarsaggr.DG_operation_hours

    # number of replacements
    replacements_number = ceil(Integer, total_DG_operation_hours/dg.lifetime) - 1
    # years that the replacements happen
    replacement_years = [i*(dg.lifetime/opervarsaggr.DG_operation_hours) for i=1:replacements_number]     # TODO verify
    # discount factors for the replacements years
    replacement_factors = [1/(1 + mg_project.discount_rate)^i for i in replacement_years]

    # component remaining life at the project end
    remaining_life = dg.lifetime - (total_DG_operation_hours - dg.lifetime * replacements_number)
    # proportional unitary salvage cost
    proportional_salvage_cost = dg.salvage_cost * remaining_life / dg.lifetime

    # present investment cost
    investment_cost = dg.investment_cost * dg.power_rated
    # present operation and maintenance cost
    om_cost = sum(dg.om_cost * dg.power_rated * opervarsaggr.DG_operation_hours * discount_factors) # depends on the nb of the DG working Hours
    # present replacement cost
    if replacements_number == 0
        replacement_cost = 0
    else
        replacement_cost = sum(dg.replacement_cost * dg.power_rated * replacement_factors)
    end
    # present salvage cost
    if remaining_life == 0
        salvage_cost = 0
    else
        salvage_cost = proportional_salvage_cost * dg.power_rated * discount_factors[mg_project.lifetime]
    end

    fuel_cost = sum(dg.fuel_cost * opervarsaggr.fuel_consumption * discount_factors)

    total_cost = investment_cost + replacement_cost + om_cost - salvage_cost + fuel_cost

    return [total_cost, investment_cost, om_cost, replacement_cost, salvage_cost, fuel_cost]
end

function annual_costs(bt::Battery, mg_project::Project, opervarsaggr::OperVarsAggr)
    if opervarsaggr.annual_throughput > 0.0
        lifetime = min(
            bt.energy_max * bt.lifetime_throughput/opervarsaggr.annual_throughput, # cycling lifetime
            bt.lifetime # calendar lifetime
        )
    else
        lifetime = bt.lifetime
    end

    c = annual_costs(
        mg_project,
        bt.energy_max,
        bt.investment_cost,
        bt.replacement_cost,
        bt.salvage_cost,
        bt.om_cost,
        0.0, 0.0,
        lifetime)
    return [c.total, c.investment, c.om, c.replacement, -c.salvage]
end

"""
    economics(mg::Microgrid, opervarsaggr::OperVarsAggr)

Return the economics results for the microgrid `mg` and
the aggregated operation variables `OperVarsAggr`.

See also: [`aggregation`](@ref)
"""
function economics(mg::Microgrid, opervarsaggr::OperVarsAggr)

    # discount factor for each year of the project
    discount_factors = [ 1/((1 + mg.project.discount_rate)^i) for i=1:mg.project.lifetime ]

    # Photovoltaic costs initialization
    PV_total_cost = 0.
    PV_investment_cost = 0.
    PV_om_cost = 0.
    PV_replacement_cost= 0.
    PV_salvage_cost = 0.
    # Wind power costs initialization
    WT_total_cost = 0.
    WT_investment_cost = 0.
    WT_om_cost = 0.
    WT_replacement_cost= 0.
    WT_salvage_cost = 0.
    # Diesel generator costs initialization
    #= DG_total_cost = 0.
    DG_investment_cost = 0.
    DG_om_cost = 0.
    DG_replacement_cost= 0.
    DG_salvage_cost = 0.
    DG_fuel_cost = 0. =#

    # NonDispatchables costs
    for i=1:length(mg.nondispatchables)
        if (typeof(mg.nondispatchables[i]) <: Photovoltaic) || (typeof(mg.nondispatchables[i]) <: PVInverter)
            PV_total_cost, PV_investment_cost, PV_om_cost, PV_replacement_cost, PV_salvage_cost = annual_costs(mg.nondispatchables[i], mg.project)
        elseif typeof(mg.nondispatchables[i]) == WindPower
            WT_total_cost, WT_investment_cost, WT_om_cost, WT_replacement_cost, WT_salvage_cost = annual_costs(mg.nondispatchables[i], mg.project)
        end
    end

    # DieselGenerator costs
    DG_total_cost, DG_investment_cost, DG_om_cost, DG_replacement_cost, DG_salvage_cost, DG_fuel_cost = annual_costs(mg.dieselgenerator, mg.project, opervarsaggr)

    # Battery costs
    BT_total_cost, BT_investment_cost, BT_om_cost, BT_replacement_cost, BT_salvage_cost = annual_costs(mg.battery, mg.project, opervarsaggr)

    # SUMMARY
    # total present investment cost
    total_investment_cost = DG_investment_cost + BT_investment_cost + PV_investment_cost + WT_investment_cost
    # total present replacement cost
    total_replacement_cost = DG_replacement_cost + BT_replacement_cost + PV_replacement_cost + WT_replacement_cost
    # total present operation and maintenance cost
    total_om_cost = DG_om_cost + BT_om_cost + PV_om_cost + WT_om_cost
    # total present salvage cost
    total_salvage_cost = DG_salvage_cost + BT_salvage_cost + PV_salvage_cost + WT_salvage_cost
    # net present cost
    npc = DG_total_cost + BT_total_cost + PV_total_cost + WT_total_cost

    # recovery factor
    recovery_factor = (mg.project.discount_rate * (1 + mg.project.discount_rate)^mg.project.lifetime)/((1 + mg.project.discount_rate)^mg.project.lifetime - 1)
    # total annualized cost
    annualized_cost = npc * recovery_factor
    # cost of energy
    coe = annualized_cost / opervarsaggr.energy_served

    # energy served over the project lifetime
    energy_served_lifetime = opervarsaggr.energy_served * sum([1.0; discount_factors[1:length(discount_factors)-1]])
    # levelized cost of energy
    lcoe = npc / energy_served_lifetime

    costs = MicrogridCosts(lcoe, coe, npc,
            total_investment_cost, total_replacement_cost, total_om_cost, total_salvage_cost,
            DG_total_cost, DG_investment_cost, DG_replacement_cost, DG_om_cost, DG_salvage_cost, DG_fuel_cost,
            BT_total_cost, BT_investment_cost, BT_replacement_cost, BT_om_cost, BT_salvage_cost,
            PV_total_cost, PV_investment_cost, PV_replacement_cost, PV_om_cost, PV_salvage_cost,
            WT_total_cost, WT_investment_cost, WT_replacement_cost, WT_om_cost, WT_salvage_cost)

    return costs
end