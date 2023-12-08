

"""
Microgrid project information

Parameters:
- lifetime (y)
- discount rate ∈ [0,1]
- time step (h)
- currency: "\$", "€"...
"""
struct Project
    "project lifetime (y)"
    lifetime::Int
    "discount rate ∈ [0,1]"
    discount_rate::Float64
    "time step (h)"
    timestep::Float64
    "currency used in price parameters and computed costs"
    currency::String
end

"""
Dispatchable power source
(e.g. Diesel generator, Gas turbine, Fuel cell with its associated tank)

All component parameters should be `Float64` except for the
*sizing parameter(s)* (here `power_rated`)
which type is parametrized as `Topt` and may be also `Float64` or
or any another `Real` type (e.g. ForwardDiff's dual number type).
"""


struct ProductionUnit{Topt<:Real}
    # Main technical parameters
    "rated power (kW)"
    power_rated::Topt
    "combustible consumption curve intercept (L/h/kW_max  ) 0for fuel cells and electrolyzer as the models used are linear " 
    consumption_intercept::Float64
    "input consumption curve slope (L/h/kW or Kg/h/kW or KW/Kg/h )"
    consumption_slope::Float64

    # Main economics parameters
    "fuel price (\$/L)"
    combustible_price::Float64
    "initial investiment price (\$/kW)"
    investment_price::Float64
    "operation & maintenance price (\$/kW/h of operation)"
    om_price_hours::Float64
    "generator lifetime (h of operation)"
    lifetime_hours::Float64

    # Secondary technical parameters (which should have a default value)
    "minimum load ratio ∈ [0,1]"
    minimum_load_ratio::Float64

    # Secondary economics parameters (which should have a default value)
    "replacement price, relative to initial investment"
    replacement_price_ratio::Float64
    "salvage price, relative to initial investment"
    salvage_price_ratio::Float64
    "fuel quantity unit (used in fuel price and consumption curve parameters)"
    input_unit::String
    output_unit::String

end
"""
Description to be added

"""

struct Tank{Topt<:Real}
    # Main technical parameters
    "rated capacity (kg or L)"
    capacity::Topt

    # Main economics parameters
    "initial investment price (\$/kg or \$/L)"
    investment_price::Float64
    "operation and maintenance price (\$/kg/y or \$/L)/y"
    om_price::Float64
    "calendar lifetime (y)"
    lifetime::Float64

    # Secondary technical parameters (which should have a default value)
    "linear loss factor α (efficiency is about 1 − α) ∈ [0,1]"
    loss_factor::Float64  # (auto - discharge rate )
    "initial level of tank ∈ [0,1]"
    ini_filling_ratio ::Float64
    "minimum level of tank ∈ [0,1]"
    min_filling_ratio::Float64
    "maximum level of tank ∈ [0,1]"
    max_filling_ratio::Float64

    # Secondary economics parameters (which should have a default value)
    "replacement price, relative to initial investment"
    replacement_price_ratio::Float64
    "salvage price, relative to initial investment"
    salvage_price_ratio::Float64
end



struct DispatchableCompound{Topt<:Real}
    generator 
    fuel_cell 

end

struct TankCompound{Topt<:Real}
    fuelTank 
    h2Tank 

end



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

struct Battery{Topt<:Real} 
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
    charge_rate::Float64 # should be type Topt in a more advanced battery model
    "max discharge power for 1 kWh (kW/kWh = h^-1)"
    discharge_rate::Float64 # should be type Topt in a more advanced battery model
    "linear loss factor α (round-trip efficiency is about 1 − 2α) ∈ [0,1]"
    loss_factor::Float64
    "minimum State of Charge ∈ [0,1]"
    SoC_min::Float64
    "maximum State of Charge ∈ [0,1]"
    SoC_max::Float64
    "initial State of Charge ∈ [0,1]"
    SoC_ini::Float64

    # Secondary economics parameters (which should have a default value)
    "replacement price, relative to initial investment"
    replacement_price_ratio::Float64
    "salvage price, relative to initial investment"
    salvage_price_ratio::Float64
end


"Base type for non-dispatchable sources (e.g. renewables like wind and solar)"
abstract type NonDispatchableSource end

"""Solar photovoltaic generator (including AC/DC converter)

See also `PVInverter` for a variant where the AC/DC converter
can be sized as well.

# About the types of the fields

All component parameters should be `Float64` except for the
*sizing parameter(s)* (here `power_rated`)
which type is parametrized as `Topt` and may be also `Float64` or
or any another `Real` type (e.g. ForwardDiff's dual number type).
"""
struct Photovoltaic{Topt<:Real} <: NonDispatchableSource
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
    derating_factor::Float64

    # Secondary economics parameters (which should have a default value)
    "replacement price, relative to initial investment"
    replacement_price_ratio::Float64
    "salvage price, relative to initial investment"
    salvage_price_ratio::Float64
end

"""Solar photovoltaic generator with adjustable AC/DC converter

See also `Photovoltaic` for a variant where the AC/DC converter
has a fixed size relative to the PV panels.

# About the types of the fields

All component parameters should be `Float64` except for the
*sizing parameter(s)* (here `power_rated`)
which type is parametrized as `Topt` and may be also `Float64` or
or any another `Real` type (e.g. ForwardDiff's dual number type).
"""
struct PVInverter{Topt<:Real} <: NonDispatchableSource
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
    derating_factor::Float64

    # Secondary economics parameters (which should have a default value)
    "replacement price, relative to initial investment"
    replacement_price_ratio::Float64
    "salvage price, relative to initial investment"
    salvage_price_ratio::Float64
end

"""Wind power generator (simple model using a given capacity factor time series)

See `capacity_from_wind` to compute capacity factor time series from wind speed.
"""
struct WindPower{Topt<:Real} <: NonDispatchableSource
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
    replacement_price_ratio::Float64
    "salvage price, relative to initial investment"
    salvage_price_ratio::Float64
end

"""Microgrid system description"""
struct Microgrid{Topt<:Real}
    "microgrid project information"
    project::Project
    "desired load time series (kW)"
    load::Vector{Float64}
    " Dispatchable Compound"
    dispatchables:: DispatchableCompound{Topt}
    "electrolyzer"
    electrolyzer
    "Tanks compound"
    tanks :: TankCompound{Topt}
    "energy storage (e.g. battery)"
    storage::Battery{Topt}
    "non-dispatchable sources (e.g. renewables like wind and solar)"
    nondispatchables
end




