# Tests for Microgrids.jl

using Microgrids
using Test

@testset "Components: PV" begin
    power_rated_PV = 5 # kW
    fPV = 0.8 # derating in [0,1]
    IT = [0., 0.5, 1.0]
    IS = 1.0
    investiment_cost_PV = 1000. # $/kW
    om_cost_PV = 20.  
    replacement_cost_PV = 1000. # $/kW
    salvage_cost_PV = 1000.
    lifetime_PV = 25
    pv = Photovoltaic(power_rated_PV, fPV, IT, IS, investiment_cost_PV, om_cost_PV, replacement_cost_PV, salvage_cost_PV, lifetime_PV)

    @test pv.power_rated == power_rated_PV
    @test production(pv) == IT * fPV * power_rated_PV;
end