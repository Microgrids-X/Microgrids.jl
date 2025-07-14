# Data types for microgrid project description, in particular each component.

"""
    SalvageType

An enum of the type of formula for salvage value calculation.

Salvage values assigns a negative cost to Microgrid components
which have a *residual lifetime* at the end of the project.

## Values

Possible values are:
- `LinearSalvage`: salvage value proportional to residual lifetime
- `ConsistentSalvage`: salvage value depends nonlinearly in the residual
  lifetime such that it is economically consistent.

Economic consistency means that, using this nonlinear formula, the NPC computation
for a given component (investment + replacement(s) - salvage) is consistent
with the annualized component cost computation (investment*CRF(component lifetime)).

Remark: zero salvage value can be obtained by setting `salvage_price_ratio=0.0`
for each Microgrid component.
"""
@enum SalvageType begin
    LinearSalvage
    ConsistentSalvage
end

"""
Microgrid project information

Parameters:
- lifetime (y)
- discount rate ∈ [0,1]
- time step (h): typically about 1.0 hour
- currency: "€" (default), "\$"...
- salvage_type: `LinearSalvage`` or `ConsistentSalvage``
"""
@kwdef mutable struct Project
    "project lifetime, i.e. economic study period (y)"
    lifetime::Int
    "discount rate ∈ [0,1]"
    discount_rate::Float64
    "time step (h)"
    timestep::Float64 = 1.0
    "currency used in price parameters and computed costs"
    currency::String = "€"
    "type of formula for salvage value calculation"
    salvage_type::SalvageType = LinearSalvage
end
# Constructors taking positional arguments with default values and conversion
Project(lifetime::Real, discount_rate::Real) = Project(
    ;lifetime=lifetime, discount_rate=discount_rate
    # default timestep, currency and salvage_type
)
Project(lifetime::Real, discount_rate::Real, timestep::Real) = Project(
    ;lifetime=lifetime, discount_rate=discount_rate, timestep=timestep
    # default currency and salvage_type
)
Project(lifetime::Real, discount_rate::Real, timestep::Real, currency::String) = Project(
    ;lifetime=lifetime, discount_rate=discount_rate, timestep=timestep, currency=currency
    # default salvage_type
)


"""
    Base.copy(p::Project)

Create a shallow copy of Microgrid project information `p`.
"""
Base.copy(p::Project) = Project(
    p.lifetime,
    p.discount_rate,
    p.timestep,
    p.currency,
    p.salvage_type
)

################################################################################
# Microgrid components
"Base type for Microgrid component (power sources, storage...)"
abstract type Component end

"test equality between two Microgrid components" # needed for *mutable* struct
Base.:(==)(c1::Comp, c2::Comp) where {Comp <: Component} = all(
    getfield(c1, field) == getfield(c2, field)
    for field in fieldnames(Comp))

"""
    Base.copy(c::Comp) where {Comp <: Component}

Create a shallow copy of a Microgrid component
"""
function Base.copy(c::Comp) where {Comp <: Component}
    T = typeof(c)
    args = (getfield(c, name) for name in fieldnames(T))
    return T(args...)
end

"""
Dispatchable power source
(e.g. Diesel generator, Gas turbine, Fuel cell)

# About the types of the fields

All component parameters should be `Float64` except for the
*sizing parameter(s)* (here `power_rated`)
which type is parametrized as `Topt` and may be also `Float64` or
or any another `Real` type (e.g. ForwardDiff's dual number type).
"""
@kwdef mutable struct DispatchableGenerator{Topt<:Real} <: Component
    # Main technical parameters
    "rated power (kW)"
    power_rated::Topt
    "fuel consumption curve intercept (L/h/kW_max)"
    fuel_intercept::Float64
    "fuel consumption curve slope (L/h/kW)"
    fuel_slope::Float64

    # Main economics parameters
    "fuel price (\$/L)"
    fuel_price::Float64
    "initial investiment price (\$/kW)"
    investment_price::Float64
    "operation & maintenance price (\$/kW/h of operation)"
    om_price_hours::Float64
    "generator lifetime (h of operation)"
    lifetime_hours::Float64

    # Secondary technical parameters (which should have a default value)
    "minimum load ratio ∈ [0,1]"
    load_ratio_min::Float64 = 0.0

    # Secondary economics parameters (which should have a default value)
    "replacement price, relative to initial investment"
    replacement_price_ratio::Float64 = 1.0
    "salvage price, relative to initial investment"
    salvage_price_ratio::Float64 = 1.0
    "fuel quantity unit (used in fuel price and consumption curve parameters)"
    fuel_unit::String = "L"
end
# Constructors taking positional arguments with default values and conversion
DispatchableGenerator(
    power_rated::Topt, fuel_intercept::Real, fuel_slope::Real, fuel_price::Real,
    investment_price::Real, om_price_hours::Real, lifetime_hours::Real
    ) where {Topt<:Real} = DispatchableGenerator(;
    power_rated, fuel_intercept, fuel_slope, fuel_price,
    investment_price, om_price_hours, lifetime_hours
)
# same but with load_ratio_min as extra argument
DispatchableGenerator(
    power_rated::Topt, fuel_intercept::Real, fuel_slope::Real, fuel_price::Real,
    investment_price::Real, om_price_hours::Real, lifetime_hours::Real,
    load_ratio_min::Real) where {Topt<:Real} = DispatchableGenerator(;
    power_rated, fuel_intercept, fuel_slope, fuel_price,
    investment_price, om_price_hours, lifetime_hours, load_ratio_min
)

"""
Battery energy storage (including AC/DC converter)

Battery dynamics is E(k+1) = E(k) − (P(k) + α|P(k)|).Δt
where α is a linear loss factor (`loss_factor` field).
It relates approximately to the roundtrip efficiency η as η = 1−2α.

# About the types of the fields

All component parameters should be `Float64` except for the
*sizing parameter(s)* (here `energy_rated`)
which type is parametrized as `Topt` and may be also `Float64` or
or any another `Real` type (e.g. ForwardDiff's dual number type).
"""
@kwdef mutable struct Battery{Topt<:Real} <: Component
    # Main technical parameters
    "rated energy capacity (kWh)"
    energy_rated::Topt

    # Main economics parameters
    "initial investment price (\$/kWh)"
    investment_price::Float64
    "operation and maintenance price (\$/kWh/y)"
    om_price::Float64
    "calendar lifetime (y)"
    lifetime_calendar::Float64
    "maximum number of cycles over life"
    lifetime_cycles::Float64

    # Secondary technical parameters (which should have a default value)
    "max charge power for 1 kWh (kW/kWh = h^-1)"
    charge_rate::Float64 = 1.0 # should be type Topt in a more advanced battery model
    "max discharge power for 1 kWh (kW/kWh = h^-1)"
    discharge_rate::Float64 = 1.0 # should be type Topt in a more advanced battery model
    "linear loss factor α (round-trip efficiency is about 1 − 2α) ∈ [0,1]"
    loss_factor::Float64 = 0.05
    "minimum State of Charge ∈ [0,1]"
    SoC_min::Float64 = 0.0
    "initial State of Charge ∈ [0,1]"
    SoC_ini::Float64 = 0.0

    # Secondary economics parameters (which should have a default value)
    "replacement price, relative to initial investment"
    replacement_price_ratio::Float64 = 1.0
    "salvage price, relative to initial investment"
    salvage_price_ratio::Float64 = 1.0
end
# Constructor taking positional arguments with default values and conversion
Battery(energy_rated::Topt, investment_price::Real, om_price::Real,
    lifetime_calendar::Real, lifetime_cycles::Real
    ) where {Topt<:Real} = Battery(;
    energy_rated, investment_price, om_price, lifetime_calendar, lifetime_cycles
)

"Base type for non-dispatchable sources (e.g. renewables like wind and solar)"
abstract type NonDispatchableSource <: Component end

"""Solar photovoltaic generator (including AC/DC converter)

See also `PVInverter` for a variant where the AC/DC converter
can be sized as well.

# About the types of the fields

All component parameters should be `Float64` except for the
*sizing parameter(s)* (here `power_rated`)
which type is parametrized as `Topt` and may be also `Float64` or
or any another `Real` type (e.g. ForwardDiff's dual number type).
"""
@kwdef mutable struct Photovoltaic{Topt<:Real} <: NonDispatchableSource
    # Main technical parameters
    "rated power (kW)"
    power_rated::Topt
    "global solar irradiance incident on the PV array time series (kW/m²)"
    irradiance::Vector{Float64}

    # Main economics parameters
    "initial investiment price (\$/kW)"
    investment_price::Float64
    "operation and maintenance price (\$/kW/y)"
    om_price::Float64
    "lifetime (y)"
    lifetime::Float64

    # Secondary technical parameters (which should have a default value)
    "derating factor (or performance ratio) ∈ [0,1]"
    derating_factor::Float64 = 0.9

    # Secondary economics parameters (which should have a default value)
    "replacement price, relative to initial investment"
    replacement_price_ratio::Float64 = 1.0
    "salvage price, relative to initial investment"
    salvage_price_ratio::Float64 = 1.0
end
# Constructors taking positional arguments with default values and conversion
Photovoltaic(power_rated::Topt, irradiance::Vector{Float64},
    investment_price::Real, om_price::Real, lifetime::Real
    ) where {Topt<:Real} =  Photovoltaic(;
    power_rated, irradiance,
    investment_price, om_price, lifetime
)
# same, with derating_factor as extra argument
Photovoltaic(power_rated::Topt, irradiance::Vector{Float64},
    investment_price::Real, om_price::Real, lifetime::Real, derating_factor::Real
    ) where {Topt<:Real} =  Photovoltaic(;
    power_rated, irradiance,
    investment_price, om_price, lifetime, derating_factor
)

"""Solar photovoltaic generator with adjustable AC/DC converter

See also `Photovoltaic` for a variant where the AC/DC converter
has a fixed size relative to the PV panels.

# About the types of the fields

All component parameters should be `Float64` except for the
*sizing parameter(s)* (here `power_rated`)
which type is parametrized as `Topt` and may be also `Float64` or
or any another `Real` type (e.g. ForwardDiff's dual number type).
"""
@kwdef mutable struct PVInverter{Topt<:Real} <: NonDispatchableSource
    # Main technical parameters
    "rated AC (inverter) power (kW)"
    power_rated::Topt
    "inverter loading ratio = ratio of DC (panels) to AC (inverter) rated powers"
    ILR::Topt
    "global solar irradiance incident on the PV array time series (kW/m²)"
    irradiance::Vector{Float64}

    # Main economics parameters
    # for AC part (inverter...)
    "initial investiment price of inverter (\$/kW_AC)"
    investment_price_ac::Float64
    "operation and maintenance price of inverter (\$/kW_AC/y)"
    om_price_ac::Float64
    "lifetime of inverter (y)"
    lifetime_ac::Float64
    # for DC part (PV panels...)
    "initial investiment price of PV panels (\$/kW_DC)"
    investment_price_dc::Float64
    "operation and maintenance price of PV panels (\$/kW_DC/y)"
    om_price_dc::Float64
    "lifetime of inverter (y)"
    lifetime_dc::Float64

    # Secondary technical parameters (which should have a default value)
    "derating factor (or performance ratio) ∈ [0,1]"
    derating_factor::Float64 = 0.9

    # Secondary economics parameters (which should have a default value)
    "replacement price, relative to initial investment"
    replacement_price_ratio::Float64 = 1.0
    "salvage price, relative to initial investment"
    salvage_price_ratio::Float64 = 1.0
end

"""Wind power generator (simple model using a given capacity factor time series)

See `capacity_from_wind` to compute capacity factor time series from wind speed.
"""
@kwdef mutable struct WindPower{Topt<:Real} <: NonDispatchableSource
    # Main technical parameters
    "rated power (kW)"
    power_rated::Topt
    "capacity factor (normalized power) time series ∈ [0,1]"
    capacity_factor::Vector{Float64}

    # Main economics parameters
    "initial investiment price (\$/kW)"
    investment_price::Float64
    "operation and maintenance price (\$/kW/y)"
    om_price::Float64
    "lifetime (y)"
    lifetime::Float64

    # Secondary economics parameters (which should have a default value)
    "replacement price, relative to initial investment"
    replacement_price_ratio::Float64 = 1.0
    "salvage price, relative to initial investment"
    salvage_price_ratio::Float64 = 1.0
end
# Constructor taking positional arguments with default values and conversion
WindPower(power_rated::Topt, capacity_factor::Vector{Float64},
    investment_price::Real, om_price::Real, lifetime::Real
    ) where {Topt<:Real} = WindPower(;
    power_rated, capacity_factor,
    investment_price, om_price, lifetime
)

"""Microgrid system description"""
@kwdef mutable struct Microgrid{Topt<:Real}
    "microgrid project information"
    project::Project
    "desired load time series (kW)"
    load::Vector{Float64}
    "dispatchable generator"
    generator::DispatchableGenerator{Topt}
    "energy storage (e.g. battery)"
    storage::Battery{Topt}
    "non-dispatchable sources (e.g. renewables like wind and solar)"
    nondispatchables #::Vector{NonDispatchableSource}
end

"""
    Base.copy(mg::Microgrid)

Create a shallow copy of Microgrid system description `mg`.

Notice: all `Microgrid` fields contains mutable data (`Project`, `Component`s
or Arrays of `Component`s). Those sub-structures are *not* copied by this method.
Copying them should be done manually afterwards if needed.

Example to make a Microgrid copy with a different discount rate:
```
mg1 = copy(mg) # assuming mg=Microgrid(...)
mg1.project = copy(mg1.project) # create a `Project` copy
mg1.project.discount_rate = 0.10 # mg.project.discount_rate unchanged
```
"""
Base.copy(mg::Microgrid) = Microgrid(
    mg.project,
    mg.load,
    mg.generator,
    mg.storage,
    mg.nondispatchables
)