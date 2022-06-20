"""
    production(photovoltaic::Photovoltaic)

Return the power output of the `photovoltaic` source.
"""
function production(photovoltaic::Photovoltaic)
    power_output = photovoltaic.derating_factor * photovoltaic.power_rated *
                   (photovoltaic.IT ./ photovoltaic.IS)
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