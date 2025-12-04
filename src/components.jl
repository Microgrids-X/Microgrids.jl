

"""
Microgrid project information
"""
struct Project
    "project lifetime (y)"
    lifetime::Int
    "discount rate ∈ [0,1]"
    discount_rate::Float32
    "time step (h)"
    timestep::Float32
    "currency used in price parameters and computed costs"
    currency::String
end

struct Sizing
   
    Cgen::Float32
    Cbatt::Float32
    Cpv::Float32
    Cwind::Float32
    Cfc::Float32
    Cel::Float32
    Htank::Float32
    Ftank::Float32
    Hb:: Float32
end

"""
ProductionUnit
Unit used to product gaz or energy
(e.g. Diesel generator, Gas turbine, Fuel cell or Electrolyzer )

All component parameters should be `Float32` except for the
*sizing parameter(s)* (here `power_rated`)
which type is parametrized as `Topt` and may be also `Float32` or
or any another `Real` type (e.g. ForwardDiff's dual number type).
Please only one O&M price should be use . Choose before a size based price (om_price) or an hourly based one om_price_hours)
"""
struct ProductionUnit{Topt<:Real}
    # Main technical parameters
    "rated power (kW)"
    power_rated::Topt
    "combustible consumption curve intercept (L/h/kW_max or kg/h/kW_max or kW_max/h/kg ); 0 for fuel cells and electrolyzer as the models used are linear " 
    consumption_intercept::Float32
    "input consumption curve slope (L/h/kW or Kg/h/kW or KW/Kg/h )"
    consumption_slope::Float32

    # Main economics parameters
    "fuel price (\$/L, \$/Kg, \$/kW)"
    combustible_price::Float32
    "initial investiment price (\$/kW)"
    investment_price::Float32
    "operation & maintenance price (\$/kW/h of operation)"
    om_price_hourly::Float32
    "operation & maintenance price (\$/kW of operation)"
    om_price::Float32
    "generator calendar lifetime "
    lifetime_calendar::Float32
    "generator lifetime (h of operation)"
    lifetime_hours::Float32
    "prod unit  maximum starts"
    lifetime_on_off :: Float32

    # Secondary technical parameters (which should have a default value)
    "minimum load ratio ∈ [0,1]"
    minimum_load_ratio::Float32

    # Secondary economics parameters (which should have a default value)
    "replacement price, relative to initial investment"
    replacement_price_ratio::Float32
    "salvage price, relative to initial investment"
    salvage_price_ratio::Float32
    "fuel quantity unit (used in fuel price and consumption curve parameters)"
    input_unit::String
    output_unit::String

end

"""
Tank( Can be a fuel or a H2 tank)
# About the types of the fields

All component parameters should be `Float32` except for the
*sizing parameter(s)* (here `capacity`)
which type is parametrized as `Topt` and may be also `Float32` or
or any another `Real` type (e.g. ForwardDiff's dual number type).
"""
struct Tank{Topt<:Real}
    # Main technical parameters
    "rated capacity (kg or L)"
    capacity::Topt

    # Main economics parameters
    "initial investment price (\$/kg or \$/L)"
    investment_price::Float32
    "operation and maintenance price (\$/kg/y or \$/L)/y"
    om_price::Float32
    "calendar lifetime (y)"
    lifetime::Float32

    # Secondary technical parameters (which should have a default value)
    "linear loss factor α (efficiency is about 1 − α) ∈ [0,1]"
    loss_factor::Float32  # 
    "initial level of tank ∈ [0,1]"
    ini_filling_ratio ::Float32
    "minimum level of tank ∈ [0,1]"
    min_filling_ratio::Float32
    "maximum level of tank ∈ [0,1]"
    max_filling_ratio::Float32

    "fuel price (\$/L, \$/Kg, \$/kW)"
    combustible_price::Float32
    
    # Secondary economics parameters (which should have a default value)
    "replacement price, relative to initial investment"
    replacement_price_ratio::Float32
    "salvage price, relative to initial investment"
    salvage_price_ratio::Float32
end

"""
Dispatchable compound
Dispatchable compound brings together the generators and fuel_cells that a `Microgrid` may contain.
"""
struct DispatchableCompound{Topt<:Real}
    "array of generators"
    generator :: Vector{ProductionUnit{Topt}}
    "array of fuel_cells"
    fuel_cell  :: Vector{ProductionUnit{Topt}}

end

"""
Tank compound
Tank compound brings together the both tanks that a `Microgrid` may contain.
We assumed that a `Microgrid` can only have one `Tank` of each type.
"""
struct TankCompound{Topt<:Real}
    "fuel tank, can be used by all diesel generators of the `Microgrid`"
    fuelTank :: Tank{Topt}
    "hydrogen tank, can be used by all Electrolyzers generators of the `Microgrid`"
    h2Tank :: Tank{Topt}

end

"""
Battery energy storage (including AC/DC converter)

Battery dynamics is E(k+1) = E(k) − (P(k) + α|P(k)|).Δt
where α is a linear loss factor (`loss_factor` field).
It relates approximately to the roundtrip efficiency η as η = 1−2α.

# About the types of the fields

All component parameters should be `Float32` except for the
*sizing parameter(s)* (here `energy_rated`)
which type is parametrized as `Topt` and may be also `Float32` or
or any another `Real` type (e.g. ForwardDiff's dual number type).
"""
struct Battery{Topt<:Real} 
    # Main technical parameters
    "rated energy capacity (kWh)"
    energy_rated::Topt
    # Main economics parameters
    "initial investment price (\$/kWh)"
    investment_price::Float32
    "operation and maintenance price (\$/kWh/y)"
    om_price::Float32
    "calendar lifetime (y)"
    lifetime_calendar::Float32
    "maximum number of cycles over life"
    lifetime_cycles::Float32

    # Secondary technical parameters (which should have a default value)
    "max charge power for 1 kWh (kW/kWh = h^-1)"
    charge_rate::Float32 # should be type Topt in a more advanced battery model
    "max discharge power for 1 kWh (kW/kWh = h^-1)"
    discharge_rate::Float32 # should be type Topt in a more advanced battery model
    "linear loss factor α (round-trip efficiency is about 1 − 2α) ∈ [0,1]"
    loss_factor::Float32
    "minimum State of Charge ∈ [0,1]"
    SoC_min::Float32
    "maximum State of Charge ∈ [0,1]"
    SoC_max::Float32
    "initial State of Charge ∈ [0,1]"
    SoC_ini::Float32

    # Secondary economics parameters (which should have a default value)
    "replacement price, relative to initial investment"
    replacement_price_ratio::Float32
    "salvage price, relative to initial investment"
    salvage_price_ratio::Float32
end


"Base type for non-dispatchable sources (e.g. renewables like wind and solar)"
abstract type NonDispatchableSource end

"""Solar photovoltaic generator (including AC/DC converter)

See also `PVInverter` for a variant where the AC/DC converter
can be sized as well.

# About the types of the fields

All component parameters should be `Float32` except for the
*sizing parameter(s)* (here `power_rated`)
which type is parametrized as `Topt` and may be also `Float32` or
or any another `Real` type (e.g. ForwardDiff's dual number type).
"""
struct Photovoltaic{Topt<:Real} <: NonDispatchableSource
    # Main technical parameters
    "rated power (kW)"
    power_rated::Topt
    "global solar irradiance incident on the PV array time series (kW/m²)"
    irradiance::Vector{Float32}

    # Main economics parameters
    "initial investiment price (\$/kW)"
    investment_price::Float32
    "operation and maintenance price (\$/kW/y)"
    om_price::Float32
    "lifetime (y)"
    lifetime::Float32

    # Secondary technical parameters (which should have a default value)
    "derating factor (or performance ratio) ∈ [0,1]"
    derating_factor::Float32

    # Secondary economics parameters (which should have a default value)
    "replacement price, relative to initial investment"
    replacement_price_ratio::Float32
    "salvage price, relative to initial investment"
    salvage_price_ratio::Float32
end

"""Solar photovoltaic generator with adjustable AC/DC converter

See also `Photovoltaic` for a variant where the AC/DC converter
has a fixed size relative to the PV panels.

# About the types of the fields

All component parameters should be `Float32` except for the
*sizing parameter(s)* (here `power_rated`)
which type is parametrized as `Topt` and may be also `Float32` or
or any another `Real` type (e.g. ForwardDiff's dual number type).
"""
struct PVInverter{Topt<:Real} <: NonDispatchableSource
    # Main technical parameters
    "rated AC (inverter) power (kW)"
    power_rated::Topt
    "inverter loading ratio = ratio of DC (panels) to AC (inverter) rated powers"
    ILR::Topt
    "global solar irradiance incident on the PV array time series (kW/m²)"
    irradiance::Vector{Float32}

    # Main economics parameters
    # for AC part (inverter...)
    "initial investiment price of inverter (\$/kW_AC)"
    investment_price_ac::Float32
    "operation and maintenance price of inverter (\$/kW_AC/y)"
    om_price_ac::Float32
    "lifetime of inverter (y)"
    lifetime_ac::Float32
    # for DC part (PV panels...)
    "initial investiment price of PV panels (\$/kW_DC)"
    investment_price_dc::Float32
    "operation and maintenance price of PV panels (\$/kW_DC/y)"
    om_price_dc::Float32
    "lifetime of inverter (y)"
    lifetime_dc::Float32

    # Secondary technical parameters (which should have a default value)
    "derating factor (or performance ratio) ∈ [0,1]"
    derating_factor::Float32

    # Secondary economics parameters (which should have a default value)
    "replacement price, relative to initial investment"
    replacement_price_ratio::Float32
    "salvage price, relative to initial investment"
    salvage_price_ratio::Float32
end

"""Wind power generator (simple model using a given capacity factor time series)

See `capacity_from_wind` to compute capacity factor time series from wind speed.
"""
struct WindPower{Topt<:Real} <: NonDispatchableSource
    # Main technical parameters
    "rated power (kW)"
    power_rated::Topt
    "capacity factor (normalized power) time series ∈ [0,1]"
    capacity_factor::Vector{Float32}

    # Main economics parameters
    "initial investiment price (\$/kW)"
    investment_price::Float32
    "operation and maintenance price (\$/kW/y)"
    om_price::Float32
    "lifetime (y)"
    lifetime::Float32

    # Secondary economics parameters (which should have a default value)
    "replacement price, relative to initial investment"
    replacement_price_ratio::Float32
    "salvage price, relative to initial investment"
    salvage_price_ratio::Float32
end

"""Microgrid system description"""
struct Microgrid{Topt<:Real}
    "microgrid project information"
    project::Project
    "desired load time series (kW)"
    load::Vector{Float32}
    " Dispatchable Compound"
    dispatchables:: DispatchableCompound{Topt}
    "electrolyzer"
    electrolyzer :: Vector{ProductionUnit{Topt}}
    "process Haber-Bosch"
    haber_bosch :: ProductionUnit{Topt}
    "Tanks compound"
    tanks :: TankCompound{Topt}
    "energy storage (e.g. battery)"
    storage::Battery{Topt}
    "non-dispatchable sources (e.g. renewables like wind and solar)"
    nondispatchables
end


