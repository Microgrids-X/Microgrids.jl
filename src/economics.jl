# Economic modeling of a microgrid project

"""Net present cost factors of some part of a Microgrid project

Cost factors can be evaluated to a *single component*
or to a *set* of components like the entire microgrid system.

Cost factors are expressed as Net Present Values, meaning they represent
*cumulated and discounted* sums over the lifetime the Microgrid project.
"""
struct CostFactors
    "total cost (initial + replacement + O&M + fuel + salvage)"
    total
    "initial investment cost"
    investment
    "replacement(s) cost"
    replacement
    "operation & maintenance (O&M) cost"
    om
    "fuel cost"
    fuel
    "salvage cost (negative)"
    salvage
end

"Add two cost factor structures, factor by factor."
function Base.:+(c1::CostFactors, c2::CostFactors)
    c = CostFactors(
        c1.total + c2.total,
        c1.investment + c2.investment,
        c1.replacement + c2.replacement,
        c1.om + c2.om,
        c1.fuel + c2.fuel,
        c1.salvage + c2.salvage
    )
    return c
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

"""
    component_costs(mg_project::Project, lifetime::Real,
        investment::Real, replacement::Real, salvage::Real,
        om_annual::Real, fuel_annual::Real)

Compute net present cost factors of a component over the Microgrid lifetime.

Cost evaluation is done from nominal and annual cost factors.

### Arguments

- mg_project: microgrid project description (e.g. with discount rate and lifetime)
- lifetime: effective lifetime of the component
- investment: initial investment cost
- replacement: nominal cost for each replacement
- salvage: nominal salvage value if component is sold at zero aging
- om_annual: nominal operation & maintenance (O&M) cost per year
- fuel_annual: nominal cost of fuel per year
"""
function component_costs(mg_project::Project, lifetime::Real,
    investment::Real, replacement::Real, salvage::Real,
    om_annual::Real, fuel_annual::Real)

    # Microgrid project parameters:
    mg_lifetime = mg_project.lifetime
    # discount factor for each year of the project
    discount_factors = [ 1/((1 + mg_project.discount_rate)^i) for i=1:mg_lifetime ]
    sum_discounts = sum(discount_factors)

    ### Operation & maintenance and fuel costs
    om_cost = om_annual * sum_discounts
    fuel_cost = fuel_annual * sum_discounts

    ### Replacement and salvage:
    if lifetime < Inf
        # number of replacements
        replacements_number = ceil(Integer, mg_lifetime/lifetime) - 1

        # net present replacement cost
        if replacements_number == 0
            replacement_cost = 0.0
        else
            # years that the replacements happen
            replacement_years = [i*lifetime for i=1:replacements_number]
            # discount factors for the replacements years
            replacement_factors = [1/(1 + mg_project.discount_rate)^i for i in replacement_years]
            replacement_cost = replacement * sum(replacement_factors)
        end

        # component remaining life at the project end
        remaining_life = lifetime*(1+replacements_number) - mg_lifetime
        # nominal *effective* salvage value (that is given remaining life)
        salvage_effective = salvage * remaining_life / lifetime

    else # Infinite lifetime (happens for components with zero usage)
        replacement_cost = 0.0
        salvage_effective = salvage # component sold "as new"
    end

    # net present salvage cost (<0)
    salvage_cost = -salvage_effective * discount_factors[mg_lifetime]

    ### Total
    total_cost = investment + replacement_cost + om_cost + fuel_cost + salvage_cost

    return CostFactors(total_cost, investment, replacement_cost, om_cost, fuel_cost, salvage_cost)
end



"""
    component_costs(nd::NonDispatchableSource, mg_project::Project)

Compute net present cost factors for a `NonDispatchableSource`.

This includes generic Photovoltaic or Wind power sources.
"""
function component_costs(nd::NonDispatchableSource, mg_project::Project)
    rating = nd.power_rated
    investment = nd.investment_price * rating
    replacement = investment * nd.replacement_price_ratio
    salvage = investment * nd.salvage_price_ratio
    om_annual = nd.om_price * rating
    fuel_annual = 0.0

    c = component_costs(
        mg_project, nd.lifetime,
        investment, replacement, salvage,
        om_annual, fuel_annual
        )
    return [c.total, c.investment, c.om, c.replacement, -c.salvage]
end

"""
    component_costs(pv::PVInverter, mg_project::Project)

Compute net present cost factors for a `PVInverter` component.
"""
function component_costs(pv::PVInverter, mg_project::Project)
    # Costs of AC part (~inverter)
    rating = pv.power_rated
    investment = pv.investment_price_ac * rating
    replacement = investment * pv.replacement_price_ratio
    salvage = investment * pv.salvage_price_ratio
    om_annual = pv.om_price_ac * rating
    fuel_annual = 0.0

    c_ac = component_costs(
        mg_project, pv.lifetime_ac,
        investment, replacement, salvage,
        om_annual, fuel_annual
        )

    # Costs of DC part (~panels)
    rating = pv.power_rated * pv.ILR # DC rated power
    investment = pv.investment_price_dc * rating
    replacement = investment * pv.replacement_price_ratio
    salvage = investment * pv.salvage_price_ratio
    om_annual = pv.om_price_dc * rating

    c_dc = component_costs(
        mg_project, pv.lifetime_dc,
        investment, replacement, salvage,
        om_annual, fuel_annual
        )
    return [c_ac.total + c_dc.total,
            c_ac.investment + c_dc.investment,
            c_ac.om + c_dc.om,
            c_ac.replacement+c_dc.replacement,
            -(c_ac.salvage+c_dc.salvage)]
end

"""
    component_costs(gen::DispatchableGenerator, mg_project::Project, oper_stats::OperationStats)

Compute net present cost factors for a `DispatchableGenerator`.
"""
function component_costs(gen::DispatchableGenerator, mg_project::Project, oper_stats::OperationStats)
    rating = gen.power_rated
    investment = gen.investment_price * rating
    replacement = investment * gen.replacement_price_ratio
    salvage = investment * gen.salvage_price_ratio
    om_annual = gen.om_price_hours * oper_stats.gen_hours * rating
    fuel_annual = gen.fuel_price * oper_stats.gen_fuel

    # effective generator lifetime in years
    lifetime = gen.lifetime_hours / oper_stats.gen_hours

    c = component_costs(
        mg_project, lifetime,
        investment, replacement, salvage,
        om_annual, fuel_annual
        )
    return [c.total, c.investment, c.om, c.replacement, -c.salvage, c.fuel]
end

"""
    component_costs(bt::Battery, mg_project::Project, oper_stats::OperationStats)

Compute net present cost factors for a `Battery`.
"""
function component_costs(bt::Battery, mg_project::Project, oper_stats::OperationStats)
    rating = bt.energy_rated
    investment = bt.investment_price * rating
    replacement = investment * bt.replacement_price_ratio
    salvage = investment * bt.salvage_price_ratio
    om_annual = bt.om_price * rating
    fuel_annual = 0.0

    if oper_stats.storage_cycles > 0.0
        lifetime = min(
            bt.lifetime_cycles/oper_stats.storage_cycles, # cycling lifetime
            bt.lifetime_calendar # calendar lifetime
        )
    else
        lifetime = bt.lifetime_calendar
    end
    c = component_costs(
        mg_project, lifetime,
        investment, replacement, salvage,
        om_annual, fuel_annual
        )

    return [c.total, c.investment, c.om, c.replacement, -c.salvage]
end

"""
    economics(mg::Microgrid, oper_stats::OperationStats)

Return the economics results for the microgrid `mg` and
the aggregated operation statistics `oper_stats`.

See also: [`aggregation`](@ref)
"""
function economics(mg::Microgrid, oper_stats::OperationStats)

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
            PV_total_cost, PV_investment_cost, PV_om_cost, PV_replacement_cost, PV_salvage_cost = component_costs(mg.nondispatchables[i], mg.project)
        elseif typeof(mg.nondispatchables[i]) == WindPower
            WT_total_cost, WT_investment_cost, WT_om_cost, WT_replacement_cost, WT_salvage_cost = component_costs(mg.nondispatchables[i], mg.project)
        end
    end

    # DieselGenerator costs
    DG_total_cost, DG_investment_cost, DG_om_cost, DG_replacement_cost, DG_salvage_cost, DG_fuel_cost = component_costs(mg.generator, mg.project, oper_stats)

    # Battery costs
    BT_total_cost, BT_investment_cost, BT_om_cost, BT_replacement_cost, BT_salvage_cost = component_costs(mg.storage, mg.project, oper_stats)

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
    coe = annualized_cost / oper_stats.served_energy

    # energy served over the project lifetime
    energy_served_lifetime = oper_stats.served_energy * sum([1.0; discount_factors[1:length(discount_factors)-1]])
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