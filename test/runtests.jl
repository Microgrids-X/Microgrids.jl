# Tests for Microgrids.jl

using Microgrids
using Test


include("economics_tests.jl")
include("optimization_tests.jl")

@testset "DispatchableGenerator" begin
    # Main parameters for DispatchableGenerator
    power_rated_gen = 10. # (kW)
    fuel_intercept = 0.05 # (L/h/kW_max)
    fuel_slope = 0.240  # (L/h/kW)
    fuel_price = 1.5 # ($/L)
    investment_price_gen = 400.
    om_price_gen = 0.02
    lifetime_gen = 15000.
    # Parameters which should have default values
    load_ratio_min = 0.33
    replacement_price_ratio = 1.2
    salvage_price_ratio = 0.8
    fuel_unit = "€"

    gen = DispatchableGenerator(power_rated_gen,
        fuel_intercept, fuel_slope, fuel_price,
        investment_price_gen, om_price_gen, lifetime_gen,
        load_ratio_min,
        replacement_price_ratio, salvage_price_ratio, fuel_unit)

    @test gen.power_rated == power_rated_gen

    @testset "DispatchableGenerator: Economics" begin
        lifetime_mg = 30.
        proj0 = Project(lifetime_mg, 0.00, 1., "€") # no discount

        # Fake operation statistics:
        oper_stats(gen_hours, gen_fuel) = OperationStats(
            0., 0., 0., 0., 0., 0.,
            1e10, gen_hours, gen_fuel,
            0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0.)

        # Case of zero usage of generator
        expected_costs_0 = [800.0, 4000.0, 0.0, 0.0,
                            3200.0, 0.0]
        component_costs(gen, proj0, oper_stats(0., 0.))
        @test round.(
            component_costs(gen, proj0, oper_stats(0., 0.));
            digits=3) == expected_costs_0
        # Case of non zero usage of generator
        expected_costs_fuel = [45812.4, 4000.0, 0.2*lifetime_mg, 0.0,
                               3200*14970/15000, 1500*lifetime_mg]
        @test round.(
            component_costs(gen, proj0, oper_stats(1., 1000.));
            digits=3 ) == expected_costs_fuel
        # TODO: add zero usage of generator of size 0...
    end
end

@testset "Photovoltaic" begin
    # Main parameters for Photovoltaic
    power_rated_pv = 10.0 # (kW)
    irradiance = [0., 0.5, 1.0] # (kW/m²)
    investment_price_pv = 1000. # ($/kW)
    om_price_pv = 20. # ($/kW/y)
    lifetime_pv = 20. # (y)
    # Parameters which should have default values
    derating_factor_pv = 0.9 # ∈ [0,1]
    replacement_price_ratio = 1.2
    salvage_price_ratio = 0.8

    pv = Photovoltaic(power_rated_pv, irradiance,
        investment_price_pv, om_price_pv,
        lifetime_pv, derating_factor_pv,
        replacement_price_ratio, salvage_price_ratio)

    @test pv.power_rated == power_rated_pv
    @test pv.lifetime == lifetime_pv
    @test production(pv) == [0.00, 0.45, 0.90] * power_rated_pv

    @testset "Photovoltaic: Economics" begin
        lifetime_mg = 30.
        proj0 = Project(lifetime_mg, 0.00, 1., "€") # no discount
        proj5 = Project(lifetime_mg, 0.05, 1., "€") # 5% discount

        expected_costs0 = [24000.00, 10000.0, 6000.00, 12000.00, 4000.00]
        expected_costs5 = [16671.65, 10000.0, 3074.49,  4522.67,  925.51]
        @test component_costs(pv, proj0) == expected_costs0
        @test round.(component_costs(pv, proj5); digits=2) == expected_costs5
    end
end

@testset "PVInverter" begin
    # Main parameters for Photovoltaic
    power_rated_pv = 5.0 # (kW AC)
    ILR = 2.0
    irradiance = [0., 0.25, 0.5, 1.0] # (kW/m²)
    # Inverter prices and lifetime:
    investment_price_ac = 100. # ($/kW)
    om_price_ac = 10.0/6 # ($/kW/y)
    lifetime_ac = 15. # (y)
    # Panel prices and lifetime:
    investment_price_dc = 1200. - investment_price_ac# ($/kW)
    om_price_dc = 20. - om_price_ac# ($/kW/y)
    lifetime_dc = 25. # (y)

    # Parameters which should have default values
    derating_factor_pv = 0.9 # ∈ [0,1]
    replacement_price_ratio = 1.0
    salvage_price_ratio = 1.0

    pvi = PVInverter(power_rated_pv, ILR, irradiance,
        investment_price_ac, om_price_ac, lifetime_ac,
        investment_price_dc, om_price_dc, lifetime_dc,
        derating_factor_pv,
        replacement_price_ratio, salvage_price_ratio)

    @test pvi.power_rated == power_rated_pv
    @test production(pvi) == [0.00, 0.45, 0.90, 1.00] * power_rated_pv

    @testset "PVInverter: Economics" begin
        lifetime_mg = 30.
        proj0 = Project(lifetime_mg, 0.00, 1., "€") # no discount
        expected_costs = [19950.0, 11500.0, 5750.0, 11500.0, 8800.0]
        @test round.(component_costs(pvi, proj0); digits=3) == expected_costs
    end
end

@testset "Economics: Battery" begin
    lifetime_mg = 25 # just a bit more than the battery
    proj0 = Project(lifetime_mg, 0.00, 1., "€") # no discount
    proj5 = Project(lifetime_mg, 0.05, 1., "€") # 5% discount

    # Main parameters for Battery
    energy_rated = 7.0 # (kWh)
    investment_price = 100.0 # ($/kWh)
    om_price = 10.0 # ($/kWh/y)
    lifetime_calendar = 20. # (y)
    lifetime_cycles = 3000.

    # Secondary parameters (which should have a default value)
    charge_rate = 2.0
    discharge_rate = 3.0
    loss_factor = 0.05
    SoC_min = 0.0
    SoC_ini = 0.5
    replacement_price_ratio = 0.9
    salvage_price_ratio = 0.8

    batt = Battery(energy_rated,
        investment_price, om_price, lifetime_calendar, lifetime_cycles,
        charge_rate, discharge_rate, loss_factor, SoC_min, SoC_ini,
        replacement_price_ratio, salvage_price_ratio)

    # Fake operation statistics for cycling:
    oper_stats(storage_cycles) = OperationStats(
        0., 0., 0., 0., 0., 0.,
        0., 0., 0.,
        storage_cycles, 0., 0., 0.,
        0., 0., 0., 0., 0., 0.)
    aggr0C   = oper_stats(0.0)
    aggr100C = oper_stats(100.0) # not enough cycles to reduce the lifetime
    aggr300C = oper_stats(300.0) # lifetime reduced to 3000/300 = 10 yr

    # no discount, with increasing amount of cycling
    @test component_costs(batt, proj0, aggr0C)   == [2660.0, 700.0, 1750.0,  630.0, 420.0]
    @test component_costs(batt, proj0, aggr100C) == [2660.0, 700.0, 1750.0,  630.0, 420.0]
    @test component_costs(batt, proj0, aggr300C) == [3430.0, 700.0, 1750.0, 1260.0, 280.0]
    # with discount
    @test round.(component_costs(batt, proj5, aggr0C); digits=2) ==  [1799.99, 700.0, 986.58, 237.44, 124.03]
end

@testset "Economics: MG" begin
    include("typical_parameters.jl")

    # Dummy time series (not used since operation simulation is bypassed)
    Pload = [1., 1., 1.]
    irradiance = [0, 0.5, 1.0]

    # Components:
    generator = DispatchableGenerator(power_rated_gen,
        fuel_intercept, fuel_slope, fuel_price,
        investment_price_gen, om_price_gen, lifetime_gen,
        load_ratio_min,
        replacement_price_ratio, salvage_price_ratio, fuel_unit)
    battery = Battery(energy_rated,
        investment_price_sto, om_price_sto, lifetime_sto, lifetime_cycles,
        charge_rate, discharge_rate, loss_factor_sto, SoC_min, SoC_ini,
        replacement_price_ratio, salvage_price_ratio)
    photovoltaic = Photovoltaic(power_rated_pv, irradiance,
        investment_price_pv, om_price_pv,
        lifetime_pv, derating_factor_pv,
        replacement_price_ratio, salvage_price_ratio)

    # Microgrid, with and without a 5% discount rate:
    proj0 = Project(lifetime, 0.0, timestep, "€")
    proj5 = Project(lifetime, discount_rate, timestep, "€")

    mg0 = Microgrid(proj0, Pload, generator, battery, [photovoltaic]);
    mg5 = Microgrid(proj5, Pload, generator, battery, [photovoltaic]);

    # Bypass of the operation simulation + aggregation:
    served_energy = 6.774979e6
    gen_hours = 3327.0
    gen_fuel = 670880.8
    storage_cycles = 209.08
    spilled_max = 4569.3
    renew_rate = 58.74029
    # Remark: 1e10 and 1e9 indicates fake values which should be nonzero
    # but do not influence economic evaluation
    oper_stats = OperationStats(
        # load stats
        served_energy, 0.0, 0.0, 0.0, 0.0, 0.0,
        # gen stats
        1e10, gen_hours, gen_fuel,
        # storage stats
        storage_cycles, 1e10, 1e10, 1e9,
        # Non-dispatchable (typ. renewables) sources stats
        1e10, spilled_max, 0.9, 1e11, 1e10, renew_rate)

    costs0 = economics(mg0, oper_stats)
    costs5 = economics(mg5, oper_stats)

    # NPC validations
    @test round(costs0.npc/1e6; digits=3) == 41.697 # M$, without discount
    @test round(costs5.npc/1e6; digits=3) == 28.353 # M$, with 5% discount
    # TODO: test LCOE
end
