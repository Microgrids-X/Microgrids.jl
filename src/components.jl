# abstract type Components end
# abstract type NonDispatchables <: Components end
abstract type NonDispatchables end

"Microgrid project information."
struct Project
    "lifetime (years)"
    lifetime::Int
    "discount rate ∈ [0,1]"
    discount_rate::Float64
    "time step (h)"
    timestep::Float64
    # TODO dispatch_type?
end

"Diesel generator parameters."
struct DieselGenerator
    "Rated power (kW)"
    power_rated::Float64   # decision variable
    "Minimum load ratio ∈ [0,1]"
    minimum_load_ratio::Float64  # ever it is on, it will work at least `min_load_ratio` of the power_max
    # min_production = min_load_ratio * power_max   # TODO - maybe it's a internal variable
    "Fuel curve intercept coefficient (L/(h × kW))"
    F0::Float64
    "Fuel curve slope (L/(h × kW))"
    F1::Float64
    
    # economics
    "Fuel cost (currency unit/L)"
    fuel_cost::Float64
    "Investiment cost (currency unit/kW)"
    investment_cost::Float64
    "Operation and maintenance cost (currency unit/(kW.h))"
    om_cost::Float64
    "Replacement cost (currency unit/kW)"
    replacement_cost::Float64
    "Salvage cost (currency unit/kW)"
    salvage_cost::Float64
    "Lifetime (h)"
    lifetime::Float64
end

"Photovoltaic parameters."
struct Photovoltaic <: NonDispatchables
    "Rated power (kW)"
    power_rated::Float64   # decision variable
    "Derating factor ∈ [0,1]"
    derating_factor::Float64
    "Incident global solar radiation (kW/m²)"
    IT::Vector{Float64}
    "Standard amount of global solar radiation (kW/m²)"
    IS::Float64

    # economics
    "Investiment cost (currency unit/kW)"
    investment_cost::Float64
    "Operation and maintenance cost (currency unit/kW)"
    om_cost::Float64
    "Replacement cost (currency unit/kW)"
    replacement_cost::Float64
    "Salvage cost (currency unit/kW)"
    salvage_cost::Float64
    "Lifetime (years)"
    lifetime::Float64

    # Photovoltaic(fPV, IT, IS, Y_PV) = new(fPV, IT, IS, Y_PV)
end

"Photovoltaic parameters with inverter issues (AC DC)"
struct PVInverter <: NonDispatchables
    "Rated power in AC (kW)"
    power_rated::Float64
    "Inverter loading ratio = PAC_rated/PDC_rated"
    ILR::Float64
    "Derating factor ∈ [0,1]"
    derating_factor::Float64
    "global solar irradiance incident on the PV array (kW/m²)"
    irradiance::Vector{Float64}

    # economics
    #AC (inverter)
    "Investiment cost of inverter (currency unit/kW)"
    investment_cost_ac::Float64
    "Operation and maintenance cost of inverter (currency unit/kW)"
    om_cost_ac::Float64
    "Replacement cost of inverter (currency unit/kW)"
    replacement_cost_ac::Float64
    "Salvage cost of inverter (currency unit/kW)"
    salvage_cost_ac::Float64
    "Lifetime of inverter (years)"
    lifetime_ac::Float64
    #DC (panels)
    "Investiment cost of pannels (currency unit/kW)"
    investment_cost_dc::Float64
    "Operation and maintenance cost of pannels (currency unit/kW)"
    om_cost_dc::Float64
    "Replacement cost of pannels (currency unit/kW)"
    replacement_cost_dc::Float64
    "Salvage cost of pannels (currency unit/kW)"
    salvage_cost_dc::Float64
    "Lifetime of pannels (years)"
    lifetime_dc::Float64
    # Photovoltaic(fPV, IT, IS, Y_PV) = new(fPV, IT, IS, Y_PV)

end

"Wind turbine parameters."
struct WindPower <: NonDispatchables
    "Rated power (kW)"
    power_rated::Float64
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
    "Investiment cost (currency unit/kW)"
    investment_cost::Float64
    "Operation and maintenance cost (currency unit/kW)"
    om_cost::Float64
    "Replacement cost (currency unit/kW)"
    replacement_cost::Float64
    "Salvage cost (currency unit/kW)"
    salvage_cost::Float64
    "Lifetime (years)"
    lifetime::Float64
end

"Battery parameters."
struct Battery
    "Initial energy (kWh)"
    energy_initial::Float64
    "Rated energy capacity (kWh)"
    energy_max::Float64    # Eb_max
    "Minimum energy level (kWh)"
    energy_min::Float64    # Eb_min  TODO - it could be the minimum state of charge too
    "Maximum charge power ∈ ``\\mathbf{R}^-`` (kW)"
    power_min::Float64     # Pb_min - charge (negative)
    "Maximum discharge power (kW)"
    power_max::Float64     # Pb_max - discharge
    "Linear loss factor ∈ [0,1]"
    loss::Float64

    # economics
    "Investiment cost (currency unit/kWh)"
    investment_cost::Float64
    "Operation and maintenance cost (currency unit/kWh)"
    om_cost::Float64
    "Replacement cost (currency unit/kWh)"
    replacement_cost::Float64
    "Salvage cost (currency unit/kWh)"
    salvage_cost::Float64
    "Lifetime (years)"
    lifetime::Float64
    "Maximum number of cycles"
    lifetime_throughput::Float64  # max throughput
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
struct Microgrid
    project::Project
    power_load::Vector{Float64}
    dieselgenerator::DieselGenerator
    # photovoltaic::Photovoltaic
    # windpower::WindPower
    battery::Battery
    nondispatchables #::Vector{NonDispatchables}
end
