# abstract type Components end
# abstract type NonDispatchables <: Components end
abstract type NonDispatchables end

"""
Microgrid project information

Parameters:
- lifetime (y)
- discount rate ∈ [0,1]
- time step (h)
- currency: "$", "€"...
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
or any another `Real` type (e.g. ForwardDiff's dual number type)
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
    "fuel price ($/L)"
    fuel_price::Float64
    "initial investiment price ($/kW)"
    investment_price::Float64
    "operation & maintenance price ($/kW/h of operation)"
    om_price_hours::Float64
    "generator lifetime (h of operation)"
    lifetime::Float64

    # Secondary technical parameters (which should have a default value)
    "minimum load ratio ∈ [0,1]"
    load_ratio_min::Float64

    # Secondary economics parameters (which should have a default value)
    "replacement price, as a fraction of initial investment price"
    replacement_price_ratio::Float64
    "salvage price, as a fraction of initial investment price"
    salvage_price_ratio::Float64
    "fuel quantity unit (used in fuel price and consumtion curve parameters)"
    fuel_unit: str = "L"
end


"""
Battery energy storage (including AC/DC converter)

Battery dynamics is E(k+1) = E(k) − (P(k) + α|P(k)|).Δt
where α is a linear loss factor (`loss` field).
It relates approximately to the roundtrip efficiency η as η = 1−2α.

# About the types of the fields

All component parameters should be `Float64` except for the
*sizing parameter(s)* (here `power_rated`)
which type is parametrized as `Topt` and may be also `Float64` or
or any another `Real` type (e.g. ForwardDiff's dual number type)
"""
struct Battery{Topt<:Real}
    "Initial energy (kWh)"
    energy_initial::Float64
    "Rated energy capacity (kWh)"
    energy_max::Topt    # Eb_max
    "Minimum energy level (kWh)"
    energy_min::Float64    # Eb_min  TODO - it could be the minimum state of charge too
    "Maximum charge power ∈ ``\\mathbf{R}^-`` (kW)"
    power_min::Topt     # Pb_min - charge (negative)
    "Maximum discharge power (kW)"
    power_max::Topt     # Pb_max - discharge
    "Linear loss factor ∈ [0,1]"
    loss::Float64

    # economics
    "initial investiment price ($/kWh)"
    investment_price::Float64
    "operation and maintenance price ($/kWh)"
    om_price::Float64
    "replacement price ($/kWh)"
    replacement_price_ratio::Float64
    "salvage price ($/kWh)"
    salvage_price_ratio::Float64
    "lifetime (y)"
    lifetime::Float64
    "maximum number of cycles"
    lifetime_throughput::Float64  # max throughput
end

"Photovoltaic parameters."
struct Photovoltaic{Topt<:Real} <: NonDispatchables
    "Rated power (kW)"
    power_rated::Topt   # decision variable
    "Derating factor ∈ [0,1]"
    derating_factor::Float64
    "Incident global solar radiation (kW/m²)"
    IT::Vector{Float64}
    "Standard amount of global solar radiation (kW/m²)"
    IS::Float64

    # economics
    "initial investiment price ($/kW)"
    investment_price::Float64
    "operation and maintenance price ($/kW)"
    om_price::Float64
    "replacement price ($/kW)"
    replacement_price_ratio::Float64
    "salvage price ($/kW)"
    salvage_price_ratio::Float64
    "lifetime (y)"
    lifetime::Float64

    # Photovoltaic(fPV, IT, IS, Y_PV) = new(fPV, IT, IS, Y_PV)
end

"Photovoltaic parameters with inverter issues (AC DC)"
struct PVInverter{Topt<:Real} <: NonDispatchables
    "Rated power in AC (kW)"
    power_rated::Topt
    "Inverter loading ratio = PAC_rated/PDC_rated"
    ILR::Topt
    "Derating factor ∈ [0,1]"
    derating_factor::Float64
    "global solar irradiance incident on the PV array (kW/m²)"
    irradiance::Vector{Float64}

    # economics
    #AC (inverter)
    "initial investiment price of inverter ($/kW)"
    investment_price_ac::Float64
    "operation and maintenance price of inverter ($/kW)"
    om_price_ac::Float64
    "replacement price of inverter ($/kW)"
    replacement_price_ratio_ac::Float64
    "salvage price of inverter ($/kW)"
    salvage_price_ratio_ac::Float64
    "lifetime of inverter (y)"
    lifetime_ac::Float64
    #DC (panels)
    "initial investiment price of pannels ($/kW)"
    investment_price_dc::Float64
    "operation and maintenance price of pannels ($/kW)"
    om_price_dc::Float64
    "replacement price of pannels ($/kW)"
    replacement_price_ratio_dc::Float64
    "salvage price of pannels ($/kW)"
    salvage_price_ratio_dc::Float64
    "lifetime of pannels (y)"
    lifetime_dc::Float64
    # Photovoltaic(fPV, IT, IS, Y_PV) = new(fPV, IT, IS, Y_PV)

end

"Wind turbine parameters."
struct WindPower{Topt<:Real} <: NonDispatchables
    "Rated power (kW)"
    power_rated::Topt
    "Cut-in speed (m/s)"
    U_cut_in::Float64
    "Cut-out speed (m/s)"
    U_cut_out::Float64
    "Rated speed (m/s)"
    U_rated::Float64
    "Wind speed at the measurement height (m/s)"
    Uanem::Vector{Float64}
    "Hub height (m)"
    zhub::Float64
    "Measurement height (m)"
    zanem::Float64
    "Roughness length (m)"
    z0::Float64
    # TODO rho
    # TODO rho0

    # economics
    "initial investiment price ($/kW)"
    investment_price::Float64
    "operation and maintenance price ($/kW)"
    om_price::Float64
    "replacement price ($/kW)"
    replacement_price_ratio::Float64
    "salvage price ($/kW)"
    salvage_price_ratio::Float64
    "lifetime (y)"
    lifetime::Float64
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

# Microgrid
struct Microgrid{Topt<:Real}
    project::Project
    power_load::Vector{Float64}
    dieselgenerator::DieselGenerator{Topt}
    # photovoltaic::Photovoltaic
    # windpower::WindPower
    battery::Battery{Topt}
    nondispatchables #::Vector{NonDispatchables}
end
