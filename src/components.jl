# abstract type Components end
# abstract type NonDispatchables <: Components end
abstract type NonDispatchables end

"Microgrid project information."
struct Project
    "lifetime (years)"
    lifetime
    "discount rate ∈ [0,1]"
    discount_rate
    "time step (h)"
    timestep
    # TODO dispatch_type?
end

"Diesel generator parameters."
struct DieselGenerator
    "Rated power (kW)"
    power_rated   # decision variable
    "Minimum load ratio ∈ [0,1]"
    minimum_load_ratio  # ever it is on, it will work at least `min_load_ratio` of the power_max
    # min_production = min_load_ratio * power_max   # TODO - maybe it's a internal variable
    "Fuel curve intercept coefficient (L/(h × kW))"
    F0
    "Fuel curve slope (L/(h × kW))"
    F1
    
    # economics
    "Fuel cost (currency unit/L)"
    fuel_cost
    "Investiment cost (currency unit/kW)"
    investment_cost
    "Operation and maintenance cost (currency unit/(kW.h))"
    om_cost
    "Replacement cost (currency unit/kW)"
    replacement_cost
    "Salvage cost (currency unit/kW)"
    salvage_cost
    "Lifetime (h)"
    lifetime
end

"Photovoltaic parameters."
struct Photovoltaic <: NonDispatchables
    "Rated power (kW)"
    power_rated   # decision variable
    "Derating factor ∈ [0,1]"
    derating_factor
    "Incident global solar radiation (kW/m²)"
    IT
    "Standard amount of global solar radiation (kW/m²)"
    IS

    # economics
    "Investiment cost (currency unit/kW)"
    investment_cost
    "Operation and maintenance cost (currency unit/kW)"
    om_cost
    "Replacement cost (currency unit/kW)"
    replacement_cost
    "Salvage cost (currency unit/kW)"
    salvage_cost
    "Lifetime (years)"
    lifetime

    # Photovoltaic(fPV, IT, IS, Y_PV) = new(fPV, IT, IS, Y_PV)
end

"Photovoltaic parameters with inverter issues (AC DC)"
struct PVInverter <: NonDispatchables
    "Rated power in AC (kW)"
    power_rated
    "Inverter loading ratio = PAC_rated/PDC_rated"
    ILR
    "Derating factor ∈ [0,1]"
    derating_factor
    "global solar irradiance incident on the PV array (kW/m²)"
    irradiance

    # economics
    #AC (inverter)
    "Investiment cost of inverter (currency unit/kW)"
    investment_cost_AC
    "Operation and maintenance cost of inverter (currency unit/kW)"
    om_cost_AC
    "Replacement cost of inverter (currency unit/kW)"
    replacement_cost_AC
    "Salvage cost of inverter (currency unit/kW)"
    salvage_cost_AC
    "Lifetime of inverter (years)"
    lifetime_AC
    #DC (panels)
    "Investiment cost of pannels (currency unit/kW)"
    investment_cost_DC
    "Operation and maintenance cost of pannels (currency unit/kW)"
    om_cost_DC
    "Replacement cost of pannels (currency unit/kW)"
    replacement_cost_DC
    "Salvage cost of pannels (currency unit/kW)"
    salvage_cost_DC
    "Lifetime of pannels (years)"
    lifetime_DC
    # Photovoltaic(fPV, IT, IS, Y_PV) = new(fPV, IT, IS, Y_PV)

end

"Wind turbine parameters."
struct WindPower <: NonDispatchables
    "Rated power (kW)"
    power_rated
    "Cut-in speed (m/s)"
    U_cut_in
    "Cut-out speed (m/s)"
    U_cut_out
    "Rated speed (m/s)"
    U_rated
    "Wind speed at the measurement height (m/s)"
    Uanem
    "Hub height (m)"
    zhub
    "Measurement height (m)"
    zanem
    "Roughness length (m)"
    z0
    # TODO rho
    # TODO rho0

    # economics
    "Investiment cost (currency unit/kW)"
    investment_cost
    "Operation and maintenance cost (currency unit/kW)"
    om_cost
    "Replacement cost (currency unit/kW)"
    replacement_cost
    "Salvage cost (currency unit/kW)"
    salvage_cost
    "Lifetime (years)"
    lifetime
end

"Battery parameters."
struct Battery
    "Initial energy (kWh)"
    energy_initial
    "Rated energy capacity (kWh)"
    energy_max    # Eb_max
    "Minimum energy level (kWh)"
    energy_min    # Eb_min  TODO - it could be the minimum state of charge too
    "Maximum charge power ∈ ``\\mathbf{R}^-`` (kW)"
    power_min     # Pb_min - charge (negative)
    "Maximum discharge power (kW)"
    power_max     # Pb_max - discharge
    "Linear loss factor ∈ [0,1]"
    loss

    # economics
    "Investiment cost (currency unit/kWh)"
    investment_cost
    "Operation and maintenance cost (currency unit/kWh)"
    om_cost
    "Replacement cost (currency unit/kWh)"
    replacement_cost
    "Salvage cost (currency unit/kWh)"
    salvage_cost
    "Lifetime (years)"
    lifetime
    "Maximum number of cycles"
    lifetime_throughput  # max throughput
end

# Operation variables - Trajectory
struct OperVarsTraj
    # load
    #= "Net load at each time instant after using the renewables power (kW)"
    Pnl_req =#
    "Net load at each time instant after dispatch (kW)"
    power_net_load
    "Unmet load/Load shedding power at each time instant (kW)"
    power_shedding
    # diesel generator
    "Diesel generator power at each time instant (kW)"
    Pgen
    # battery
    "Battery energy at each time instant (kWh)"
    Ebatt
    "Battery power at each time instant (kW)"
    Pbatt
    "Maximum battery discharge power at time t (kW)"
    Pbatt_dmax
    "Maximum battery charge power at time t (kW)"
    Pbatt_cmax
    # renewables sources
    "Renewables curtailment power at each time instant (kW)"
    power_curtailment
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
    power_load
    dieselgenerator::DieselGenerator
    # photovoltaic::Photovoltaic
    # windpower::WindPower
    battery::Battery
    nondispatchables::Vector{NonDispatchables}
end
