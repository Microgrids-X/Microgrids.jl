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

@testset "Economics: PV" begin
    lifetime_mg = 30
    discount_rate = 0.05

    power_rated_PV = 10 # kW
    fPV = 1.
    IT = [0., 0.5, 1.0]
    IS = 1.0
    investiment_cost_PV = 1000.
    om_cost_PV = 20.  
    replacement_cost_PV = 1200.
    salvage_cost_PV = 800.
    lifetime_PV = 20

    project0 = Project(lifetime_mg, 0.00, 1.) # no discount
    project5 = Project(lifetime_mg, 0.05, 1.) # 5% discount

    pv = Photovoltaic(power_rated_PV, fPV, IT, IS, investiment_cost_PV, om_cost_PV, replacement_cost_PV, salvage_cost_PV, lifetime_PV)

    @test annual_costs(pv, project0) ==                   [24000.0,  10000.0, 6000.0, 12000.0, 4000.0]
    @test round.(annual_costs(pv, project5); digits=2) == [16671.65, 10000.0, 3074.49, 4522.67, 925.51]
end