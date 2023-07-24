"""
    production(photovoltaic::Photovoltaic)

Return the power output of the `photovoltaic` source.
"""
function production(photovoltaic::Photovoltaic{T})::Vector{T} where {T<:Real}
    power_output = photovoltaic.derating_factor * photovoltaic.power_rated *
                   photovoltaic.irradiance
    return power_output
end

"""
    production(photovoltaic::Photovoltaic, k::Int)

Return the power output of the `photovoltaic` source at instant `k`.
"""
function production(photovoltaic::Photovoltaic{T}, k::Int)::T where {T<:Real}
    power_output = photovoltaic.derating_factor * photovoltaic.power_rated *
                   photovoltaic.irradiance[k]
    return power_output
end

"""
    production(photovoltaic::PVInverter)

Return the power output of the `photovoltaic` source.
"""
function production(photovoltaic::PVInverter{T})::Vector{T} where {T<:Real}
    PDC_rated=photovoltaic.ILR*photovoltaic.power_rated
    PAC=zeros(T, length(photovoltaic.irradiance))
    for i=1:length(photovoltaic.irradiance)
        PDC_i=photovoltaic.irradiance[i]*PDC_rated*photovoltaic.derating_factor
        PAC[i]=min(PDC_i,photovoltaic.power_rated)
    end
    return PAC
end

"""
    production(photovoltaic::PVInverter, k::Int)

Return the power output of the `photovoltaic` source at instant `k`.
"""
function production(photovoltaic::PVInverter{T}, k::Int)::T where {T<:Real}
    PDC_rated=photovoltaic.ILR*photovoltaic.power_rated
    PDC=photovoltaic.irradiance[k]*PDC_rated*photovoltaic.derating_factor
    PAC=min(PDC,photovoltaic.power_rated)
    return PAC
end

"""
    production(wind::WindPower)

Return the power output of the `wind` source. A quadratic approximation is used
for the power curve.

!!! note

    It is defined only for the non-dispatchable sources.
"""
function production(windpower::WindPower)
    Uhub = windpower.Uanem * log(windpower.zhub/windpower.z0) / log(windpower.zanem/windpower.z0)
    power_output = zeros(Real,length(windpower.Uanem))   # choisir le type que le ForwardDiff necessite
    for i=1:length(windpower.Uanem)
        Pwind_aux = 0
        if windpower.U_cut_in < Uhub[i] < windpower.U_rated
            Pwind_aux = windpower.power_rated * ((Uhub[i] - windpower.U_cut_in)/(windpower.U_rated - windpower.U_cut_in))^2  # usando o modelo quadrático
        elseif windpower.U_rated <= Uhub[i] <= windpower.U_cut_out
            Pwind_aux = windpower.power_rated
        end
        power_output[i] = Pwind_aux
    end
    return power_output
end