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
    proj0 = Project(lifetime_mg, 0.00, 1.) # no discount
    proj5 = Project(lifetime_mg, 0.05, 1.) # 5% discount

    power_rated_PV = 10 # kW
    fPV = 1.
    IT = [0., 0.5, 1.0]
    IS = 1.0
    investiment_cost_PV = 1000.
    om_cost_PV = 20.  
    replacement_cost_PV = 1200.
    salvage_cost_PV = 800.
    lifetime_PV = 20

    pv = Photovoltaic(power_rated_PV, fPV, IT, IS, investiment_cost_PV, om_cost_PV, replacement_cost_PV, salvage_cost_PV, lifetime_PV)

    @test annual_costs(pv, proj0) ==                   [24000.0,  10000.0, 6000.0, 12000.0, 4000.0]
    @test round.(annual_costs(pv, proj5); digits=2) == [16671.65, 10000.0, 3074.49, 4522.67, 925.51]
end

@testset "Economics: Battery" begin
    lifetime_mg = 25 # just a bit more than the battery
    proj0 = Project(lifetime_mg, 0.00, 1.) # no discount
    proj5 = Project(lifetime_mg, 0.05, 1.) # 5% discount
    
    energy_initial = 0.
    energy_max = 7
    energy_min = 0
    power_min = -1.0*energy_max
    power_max = +1.0*energy_max
    loss = 0.05
    investiment_cost_BT = 100.
    om_cost_BT = 10
    replacement_cost_BT = 90.
    salvage_cost_BT = 80
    lifetime_BT = 20
    lifetime_thrpt = 3000

    batt = Battery(energy_initial, energy_max, energy_min, power_min, power_max, loss, investiment_cost_BT, om_cost_BT, replacement_cost_BT, salvage_cost_BT, lifetime_BT, lifetime_thrpt)

    # Fake operation data: 1k cycles/year, 
    aggr0C   = OperVarsAggr(0., 0., 0., 0., 0., 0., 0.,    0., 0., 0.)
    aggr100C = OperVarsAggr(0., 0., 0., 0., 0., 0., 0., 100.0*energy_max, 0., 0.) # not enough cycles to reduce the lifetime
    aggr300C = OperVarsAggr(0., 0., 0., 0., 0., 0., 0., 300.0*energy_max, 0., 0.) # lifetime reduced to 3000/300 = 10 yr

    # no discount, with increasing amount of cycling
    @test annual_costs(batt, proj0, aggr0C)   == [2660.0, 700.0, 1750.0,  630.0, 420.0]
    @test annual_costs(batt, proj0, aggr100C) == [2660.0, 700.0, 1750.0,  630.0, 420.0]
    @test annual_costs(batt, proj0, aggr300C) == [3430.0, 700.0, 1750.0, 1260.0, 280.0]
    # with discount
    @test round.(annual_costs(batt, proj5, aggr0C); digits=2) ==  [1799.99, 700.0, 986.58, 237.44, 124.03]
end

@testset "Economics: MG" begin
    Pload = [1., 1., 1.]
    
    lifetime = 25 # yr
    discount_rate = 0.05
    timestep = 1 # h

    proj0 = Project(lifetime, 0, timestep)
    proj5 = Project(lifetime, discount_rate, timestep)

    power_rated_PV = 6000 # kW
    fPV = 1.
    IT = [0, 0.5, 1.0]
    IS = 1.0
    investiment_cost_PV = 1200.
    om_cost_PV = 20.  
    replacement_cost_PV = 1200.
    salvage_cost_PV = 1200.
    lifetime_PV = 25

    photovoltaic = Photovoltaic(power_rated_PV, fPV, IT, IS, investiment_cost_PV, om_cost_PV, replacement_cost_PV, salvage_cost_PV, lifetime_PV)

    energy_initial = 0.
    energy_max = 9000
    energy_min = 0
    power_min = -1.0*energy_max
    power_max = +1.0*energy_max
    loss = 0.05
    investiment_cost_BT = 350.
    om_cost_BT = 10.
    replacement_cost_BT = 350.
    salvage_cost_BT = 350.
    lifetime_BT = 15
    lifetime_thrpt = 3000

    battery = Battery(energy_initial, energy_max, energy_min, power_min, power_max, loss, investiment_cost_BT, om_cost_BT, replacement_cost_BT, salvage_cost_BT, lifetime_BT, lifetime_thrpt)

    power_rated_DG = 1800.
    min_load_ratio = 0
    F0 = 0.0
    F1 = 0.240
    fuel_cost = 1.
    investiment_cost_DG = 400.
    om_cost_DG = 0.02
    replacement_cost_DG = 400.
    salvage_cost_DG = 400.
    lifetime_DG = 15000

    dieselgenerator = DieselGenerator(power_rated_DG, min_load_ratio, F0, F1, fuel_cost, investiment_cost_DG, om_cost_DG, replacement_cost_DG, salvage_cost_DG, lifetime_DG)

    mg0 = Microgrid(proj0, Pload, dieselgenerator, battery, [photovoltaic]);
    mg5 = Microgrid(proj5, Pload, dieselgenerator, battery, [photovoltaic]);

    # Bypass of the operation simulation + aggregation:
    opervarsaggr = OperVarsAggr(6.774979e6, 0, 0, 0, 0.0, 3327, 670880.8054857135, 1.881753568922309e6, 4569.3, 58.74028997693115)

    costs0 = economics(mg0, opervarsaggr)
    costs5 = economics(mg5, opervarsaggr)
    # NPC validations
    @test round(costs0.npc/1e6; digits=3) == 41.697 # M$, without discount
    @test round(costs5.npc/1e6; digits=3) == 28.353 # M$, with 5% discount
    # TODO: test LCOE

end
