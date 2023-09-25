"""
    production(photovoltaic::Photovoltaic)

Return the power output of the `photovoltaic` source.
"""
function production(photovoltaic::Photovoltaic)
    power_output = photovoltaic.derating_factor * photovoltaic.power_rated *
                   photovoltaic.irradiance
    return power_output
end

"""
    production(photovoltaic::PVInverter)

Return the power output of the `photovoltaic` source.
"""
function production(photovoltaic::PVInverter)
    PDC_rated=photovoltaic.ILR*photovoltaic.power_rated
    T = typeof(PDC_rated)
    PAC=zeros(T, length(photovoltaic.irradiance))
    for i=1:length(photovoltaic.irradiance)
        PDC_i=photovoltaic.irradiance[i]*PDC_rated*photovoltaic.derating_factor
        PAC[i]=min(PDC_i,photovoltaic.power_rated)
    end
    return PAC
end

"""
    capacity_from_wind(v::Real; TSP::Real, Cp=0.50, v_out=25.0, α=3.0)

Compute capacity factor (normalized power) of a wind turbine,
using a generic parametrized power curve P(v), for a given a wind speed `v` (m/s).

Model parameters are:
- Turbine Specific Power `TSP`, in W/m², typically 200 – 400.
- Maximum power coefficient `Cp` (used before saturation)
  should be smaller than Betz' limit of 16/27.
- Cut-out wind speed is `v_out` (m/s).

A fixed Cp model is used, with a soft saturation when reaching maximum power.
This soft saturation, based on LogSumExp, is tuned with `α`
(default: 3.0, higher yields sharper transition).

Air is assumed to have fixed density ρ=1.225 kg/m³.
"""
function capacity_from_wind(v::Real; TSP::Real, Cp=0.50, v_out=25.0, α=3.0)
    ρ = 1.225 # kg/m³ at 15°C
    # Normalized power from the wind, without saturation:
    cf = 0.5*Cp*ρ/TSP * v^3
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

"""
    production(windpower::WindPower)

Return the power output of the `windpower` source.
"""
function production(windpower::WindPower)
    power_output = windpower.power_rated * windpower.capacity_factor
    return power_output
end