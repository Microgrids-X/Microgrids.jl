

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
(e.g. Diesel generator, Gas turbine, Fuel cell)

# About the types of the fields

All component parameters should be `Float64` except for the
*sizing parameter(s)* (here `power_rated`)
which type is parametrized as `Topt` and may be also `Float64` or
or any another `Real` type (e.g. ForwardDiff's dual number type).
"""
struct DispatchableGenerator{Topt<:Real}
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
    load_ratio_min::Float64

    # Secondary economics parameters (which should have a default value)
    "replacement price, relative to initial investment"
    replacement_price_ratio::Float64
    "salvage price, relative to initial investment"
    salvage_price_ratio::Float64
    "fuel quantity unit (used in fuel price and consumption curve parameters)"
    fuel_unit::String
end


"""
Battery energy storage (including AC/DC converter)

Battery dynamics is E(k+1) = E(k) − (P(k) + α|P(k)|).Δt
where α is a linear loss factor (`loss_factor` field).
It relates approximately to the roundtrip efficiency η as η = 1−2α.

# About the types of the fields

All component parameters should be `Float64` except for the
*sizing parameter(s)* (here `energy_rated` and optionally `(dis)charge_rate`)
which type is parametrized as `Topt1/2` and may be also `Float64` or
or any another `Real` type (e.g. ForwardDiff's dual number type).
"""
struct Battery{Topt1<:Real,Topt2<:Real}
    # Main technical parameters
    "rated energy capacity (kWh)"
    energy_rated::Topt1

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
    charge_rate::Topt2
    "max discharge power for 1 kWh (kW/kWh = h^-1)"
    discharge_rate::Topt2
    "linear loss factor α (round-trip efficiency is about 1 − 2α) ∈ [0,1]"
    loss_factor::Float64
    "minimum State of Charge ∈ [0,1]"
    SoC_min::Float64
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

"""Solar photovoltaic generator with ajdustable AC/DC converter

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

function capacity_from_wind(v, TSP, Cp=0.50, v_out=25.0, α=3.0)
    """Compute capacity factor (normalized power) of a wind turbine,
    using a generic parametrized power curve P(v),
    for a given a wind speed `v` (m/s).

    Model parameters are:
    Turbine Specific Power `TSP`, in W/m², typically 200 – 400.
    Maximum power coefficient `Cp` (used before saturation)
    should be smaller than Betz' limit of 16/27.
    Cut-out wind speed is `v_out` (m/s).

    A fixed Cp model is used, with a soft saturation when reaching maximum power.
    This soft saturation, based on LogSumExp, is tuned with `α`
    (default: 3.0, higher yields sharper transition).

    Air is assumed to have fixed density ρ=1.225 kg/m³.
    """
    ρ = 1.225 # kg/m³ at 15°C
    # Normalized power from the wind, without saturation:
    cf = 0.5*Cp*ρ/TSP * v.^3
    # saturation using a smooth min based on LogSumExp
    cf = -log(exp(-α) + exp(-α*cf)) / α
    # saturate negative values (due to the smooth min)
    if cf < 0.0
        cf = 0.0
    end
    # Cut-out wind speed:
    if v > v_out
        cf = 0.0
    end
    return cf
end

# Operation variables - Trajectory
struct OperVarsTraj{T<:Real}
    # load
    #= "Net load at each time instant after using the renewables power (kW)"
    Pnl_req =#
    "Net load at each time instant after dispatch (kW)"
    power_net_load::Vector{T}
    "Unmet load/Load shedding power at each time instant (kW)"
    power_shedding::Vector{T}
    # diesel generator
    "Diesel generator power at each time instant (kW)"
    Pgen::Vector{T}
    # battery
    "Battery energy at each time instant (kWh)"
    Ebatt::Vector{T}
    "Battery power at each time instant (kW)"
    Pbatt::Vector{T}
    "Maximum battery discharge power at time t (kW)"
    Pbatt_dmax::Vector{T}
    "Maximum battery charge power at time t (kW)"
    Pbatt_cmax::Vector{T}
    # renewables sources
    "Renewables curtailment power at each time instant (kW)"
    power_curtailment::Vector{T}
end

# Operation variables - Aggregation
struct OperVarsAggr
    # load
    "Load energy served in one year (kWh)"
    energy_served
    "Maximum load shedding power (kW)"
    power_shedding_max
    "Maximum consecutive duration of load shedding (h)"
    shedding_duration_max
    "Load shedding energy in one year (kWh)"
    energy_shedding_total
    "Ratio between load shedding energy and total consumption in one year (%)"
    shedding_rate
    # diesel generator
    "Number of diesel generator operation hours in one year (h)"
    DG_operation_hours
    "Diesel consumption in one year (L)"
    fuel_consumption
    # battery
    "Number of completed battery cycles in one year"
    annual_throughput
    # renewables sources
    "Maximum renewables curtailment power (kW)"
    power_curtailment_max
    "Ratio between energy supplied by renewables and energy served in one year (%)"
    renewables_rate
end

"""Microgrid system description"""
struct Microgrid{Topt<:Real}
    "microgrid project information"
    project::Project
    "desired load time series (kW)"
    load::Vector{Float64}
    "dispatchable generator"
    generator::DispatchableGenerator{Topt}
    "energy storage (e.g. battery)"
    storage::Battery{Topt,Topt}
    "non-dispatchable sources (e.g. renewables like wind and solar)"
    nondispatchables #::Vector{NonDispatchableSource}
end
