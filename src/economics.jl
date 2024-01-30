# Economic modeling of a microgrid project

### Structures to hold costs

"""Net present cost factors of some part of a Microgrid project

Cost factors can be evaluated to a *single component*
or to a *set* of components like the entire microgrid system.

Cost factors are expressed as Net Present Values, meaning they represent
*cumulated and discounted* sums over the lifetime the Microgrid project.
"""
@kwdef struct CostFactors
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

# Arithmetic for CostFactors: +,*,/ and round

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

"Multiply cost factor `c` by scalar `a`"
function Base.:*(c::CostFactors, a::Real)
    c = CostFactors(
        c.total * a,
        c.investment * a,
        c.replacement * a,
        c.om * a,
        c.fuel * a,
        c.salvage * a
    )
    return c
end

"Multiply scalar `a` by cost factor `c`"
Base.:*(a::Real, c::CostFactors) = *(c, a) # commutativity

"Divide cost factor `c` by scalar `a`"
Base.:/(c::CostFactors, a::Real) = *(c, inv(a)) # commutativity

"""
    round(c::CostFactors, [r::RoundingMode])
    round(c::CostFactors, [r::RoundingMode]; digits::Integer=0)
    round(c::CostFactors, [r::RoundingMode]; sigdigits::Integer)

Round each field of `CostFactors` object `c`.

Notice that rounding is done on each field separately so that
the rounded `total` field may end up being different from the sum
of all the other fields.
"""
function Base.round(c::CostFactors, r=RoundNearest;
                    digits::Union{Nothing,Integer}=nothing,
                    sigdigits::Union{Nothing,Integer}=nothing)
    if sigdigits === nothing # round on digits
        if digits === nothing
            digits=0
        end
        cr = CostFactors(
            round(c.total, r; digits=digits),
            round(c.investment, r; digits=digits),
            round(c.replacement, r; digits=digits),
            round(c.om, r; digits=digits),
            round(c.fuel, r; digits=digits),
            round(c.salvage, r; digits=digits)
        )
    else # round on sigdigits
        if digits !== nothing
            throw(ArgumentError("`round` cannot use both `digits` and `sigdigits` arguments."))
        end
        cr = CostFactors(
            round(c.total, r; sigdigits=sigdigits),
            round(c.investment, r; sigdigits=sigdigits),
            round(c.replacement, r; sigdigits=sigdigits),
            round(c.om, r; sigdigits=sigdigits),
            round(c.fuel, r; sigdigits=sigdigits),
            round(c.salvage, r; sigdigits=sigdigits)
        )
    end
    return cr
end


"""Cost factors of each component of a Microgrid

also includes `system` cost (all components) and two key economic data:
`npc` and `lcoe`:
- `npc` is equal to `system.total`
- `lcoe` is `npc` divided by the discounted sum of served energy,
  or, assuming a constant annual served energy, it is also
  annualized npc divided by yearly served energy
"""
@kwdef struct MicrogridCosts
    "levelized cost of electricity (\$/kWh)"
    lcoe:: Real
    "net present cost of the microgrid (\$)"
    npc:: Real
    "costs of all components"
    system:: CostFactors
    "costs of generator"
    generator:: CostFactors
   
    "costs of energy storage"
    storage:: CostFactors
    "costs of each non-dispatchable source"
    nondispatchables:: Vector{CostFactors}
    electrolyzer:: CostFactors
    fuel_cell:: CostFactors
    fuel_tank:: CostFactors
    h2_tank:: CostFactors
end

### component_costs methods

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
    component_costs(gen::ProductionUnit, mg_project::Project, oper_stats::OperationStats)

Compute net present cost factors for a `ProductionUnit`.
"""
function component_costs(prod_unit::ProductionUnit, mg_project::Project, prod_unit_hours, prod_unit_cons)
    rating = prod_unit.power_rated
    investment = prod_unit.investment_price * rating
    replacement = investment * prod_unit.replacement_price_ratio
    salvage = investment * prod_unit.salvage_price_ratio
    om_annual = prod_unit.om_price_hours * rating
    fuel_annual = prod_unit.combustible_price * prod_unit_cons

    # effective generator lifetime (in years)
    lifetime = prod_unit.lifetime_hours / prod_unit_hours

    c = component_costs(
        mg_project, lifetime,
        investment, replacement, salvage,
        om_annual, fuel_annual
        )
    return c
end

"""
    component_costs(tank::Tank, mg_project::Project)

Compute net present cost factors for a `Tank`.
"""
function component_costs(tank::Tank, mg_project::Project)
    rating = tank.capacity
    investment = tank.investment_price * rating + (tank.ini_filling_ratio * tank.combustible_price * rating )
    replacement = investment * tank.replacement_price_ratio
    salvage = investment * tank.salvage_price_ratio
    om_annual = tank.om_price * rating
    fuel_annual = 0.0

    c = component_costs(
        mg_project, tank.lifetime,
        investment, replacement, salvage,
        om_annual, fuel_annual
        )
    return c
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

    # effective battery lifetime (in years)
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

    return c
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
    return c
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

    c = c_ac + c_dc
    return c
end

### Economic evaluation of an entire microgrid project

"""
    economics(mg::Microgrid, oper_stats::OperationStats)

Return the economics results for the microgrid `mg` and
the aggregated operation statistics `oper_stats`.

See also: [`aggregation`](@ref)
"""
function economics(mg::Microgrid, oper_stats::OperationStats)
    # Dispatchable generator
    gen_costs = component_costs(mg.dispatchables.generator[1], mg.project, oper_stats.gen_hours,oper_stats.gen_fuel)
    elyz_costs = component_costs(mg.electrolyzer[1], mg.project, oper_stats.elyz_hours,oper_stats.elyz_consumed_energy)
    fc_costs = component_costs(mg.dispatchables.fuel_cell[1], mg.project, oper_stats.fc_hours,oper_stats.h2_consumed)
    fuel_tantk_cost = component_costs(mg.tanks.fuelTank, mg.project)
    h2_tantk_cost = component_costs(mg.tanks.h2Tank, mg.project)
    # Energy storage
    sto_costs = component_costs(mg.storage, mg.project, oper_stats)

    # Non-dispatchable sources (e.g. renewables like wind and solar)

    # NonDispatchables costs
    nd_costs = collect(
        component_costs(mg.nondispatchables[i], mg.project)
        for i in 1:length(mg.nondispatchables)
        )

    # Capital recovery factor (CRF)
    discount_factors = [1/((1 + mg.project.discount_rate)^i)
                        for i = 1:mg.project.lifetime]
    crf = 1/sum(discount_factors)
    # Cost of all components and NPC of the project
    system_costs = gen_costs + sto_costs + sum(nd_costs) + elyz_costs + fc_costs + fuel_tantk_cost + h2_tantk_cost
    npc = system_costs.total
    # levelized cost of energy
    annualized_cost = npc*crf # $/y
    lcoe = annualized_cost / oper_stats.served_energy # ($/y) / (kWh/y) → $/kWh

    costs = MicrogridCosts(lcoe, npc,
        system_costs, gen_costs, sto_costs, nd_costs, elyz_costs,fc_costs,fuel_tantk_cost,h2_tantk_cost
    )

    return costs
end