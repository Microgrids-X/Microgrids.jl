# Tests for Microgrid components creation and their cost calculation

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
        c_0 = CostFactors(
            800.0, 4000.0, 0.0, 0.0, 0.0, -3200.0
        )
        @test round(
            component_costs(gen, proj0, oper_stats(0., 0.));
            digits=3) == c_0
        # Case of non zero usage of generator
        c_fuel = CostFactors(
            total=45812.4,
            investment=4000.0,
            replacement=0.0,
            om=0.2*lifetime_mg,
            fuel=1500*lifetime_mg,
            salvage=-3200*14970/15000
        )
        @test round(
            component_costs(gen, proj0, oper_stats(1., 1000.));
            digits=3 ) == c_fuel
        # TODO: add zero usage of generator of size 0...
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

    # Costs with no discount, with increasing amount of cycling
    c_0C = CostFactors( # cost with no to little cycling (calendar lifetime dominant)
        total=2660.0,
        investment=700.0,
        replacement=700.0*0.9, # 630.0
        om=70.0*25, # 1750.0
        fuel=0.0,
        salvage=-420.0
    )
    c_300C = CostFactors( # cost with lifetime dominated by cycling
        total=3430.0,
        investment=700.0,
        replacement=700.0*0.9*2, # 1260.0
        om=70.0*25, # 1750.0
        fuel=0.0,
        salvage=-280.0
    )
    @test component_costs(batt, proj0, aggr0C)   == c_0C
    @test component_costs(batt, proj0, aggr100C) == c_0C # 100 c/y is like 0
    @test component_costs(batt, proj0, aggr300C) == c_300C
    # Costs *with* discount
    c5_0C = CostFactors(
        total=1799.99,
        investment=700.0,
        replacement=237.44,
        om=986.58,
        fuel=0.0,
        salvage=-124.03
    )
    @test round(component_costs(batt, proj5, aggr0C); digits=2) == c5_0C
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

        c0 = CostFactors(
            total=24000.00,
            investment=10000.0,
            replacement=12000.00,
            om=6000.00,
            fuel=0.0,
            salvage=-4000.00
        )
        c5 = CostFactors(
            total=16671.65,
            investment=10000.0,
            replacement=4522.67,
            om=3074.49,
            fuel=0.0,
            salvage=-925.51
        )
        @test component_costs(pv, proj0) == c0
        @test round(component_costs(pv, proj5); digits=2) == c5
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
        c = CostFactors(
            total=19950.0,
            investment=11500.0,
            replacement=11500.0,
            om=5750.0,
            fuel=0.0,
            salvage=-8800.0
        )
        @test round(component_costs(pvi, proj0); digits=3) == c
    end
end

@testset "WindPower" begin
    # Wind turbine parameters fitted to an EWT 900 kW DW52
    S_D52 = pi * (52/2)^2 # rotor swept area m²
    TSP_D52 = 900e3/S_D52 # W/m²
    v_out = 25.0 # cut-out speed m/s
    Cp_D52, α_D52 = 0.521, 3.1; # fitted from actual power curve

    wind_speed = [0., 2., 3.,    5.,    7.,    10.,   15., 25., 25.1] # m/s
    cf_exp =     [0., 0., 0.005, 0.075, 0.227, 0.630, 0.997, 1.0, 0.]
    cf_wind = capacity_from_wind.(wind_speed; TSP=TSP_D52, Cp=Cp_D52, v_out=v_out, α=α_D52)

    @test cf_wind ≈ cf_exp atol=1e-3

    # Main parameters for WindPower
    power_rated_wind = 1000. # rated power (kW)
    investment_price_wind = 3000. # initial investiment price ($/kW)
    om_price_wind = 60.# operation and maintenance price ($/kW/y)
    lifetime_wind = 25. # lifetime (y)

    # Parameters which should have default values
    replacement_price_ratio = 1.0
    salvage_price_ratio = 1.0

    windgen = WindPower(power_rated_wind, cf_wind,
        investment_price_wind, om_price_wind, lifetime_wind,
        replacement_price_ratio, salvage_price_ratio)

    # Wind power production time series
    @test production(windgen) ≈ cf_exp*power_rated_wind atol=1

    @testset "WindPower: Economics" begin
        lifetime_mg = 30.
        proj0 = Project(lifetime_mg, 0.00, 1., "€") # no discount

        c0 = CostFactors(
            total=5.4e6,
            investment=3.0e6,
            replacement=3.0e6,
            om=60e3*lifetime_mg,
            fuel=0.0,
            salvage=-3.0e6*20/25
        )
        @test round(component_costs(windgen, proj0); digits=3) == c0
    end
end