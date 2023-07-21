# Tests for economics calculation
# Remark: some test for the economics of components are grouped along
# the tests of their respective component.

@testset "Economics: Cost factors" begin
    investment1 = 1.0
    replacement1 = 2.0
    om1 = 3.0
    fuel1 = 4.0
    salvage1 = -0.5
    total1 = investment1 + replacement1 + om1 + fuel1 + salvage1
    c1 =  CostFactors(total1, investment1, replacement1, om1, fuel1, salvage1)

    investment2 = 10.0
    replacement2 = 20.0
    om2 = 30.0
    fuel2 = 40.0
    salvage2 = -5.0
    total2 = investment2 + replacement2 + om2 + fuel2 + salvage2
    c2 =  CostFactors(total2, investment2, replacement2, om2, fuel2, salvage2)

    # Addition
    ctot = c1 + c2
    @test ctot.total == c1.total + c2.total
    @test ctot.investment == c1.investment + c2.investment
    @test ctot.replacement == c1.replacement + c2.replacement
    @test ctot.om == c1.om + c2.om
    @test ctot.fuel == c1.fuel + c2.fuel
    @test ctot.salvage == c1.salvage + c2.salvage

    # sum using the sum function (useful to sum costs of nondispatchables)
    @test sum([c1, c2]) == ctot
    @test sum([c1, c2, c2]) != ctot

    # Multiplication with scalar
    @test c1+c1 == 2*c1
    @test c1+c1 == c1*2
    # Division with scalar
    @test (c1+c1)/2 == c1

    # Cost rounding
    investment3 = 1.11
    replacement3 = 2.22
    om3 = 3.33
    fuel3 = 4.44
    salvage3 = -0.555
    total3 = investment3 + replacement3 + om3 + fuel3 + salvage3
    c3 =  CostFactors(total3, investment3, replacement3, om3, fuel3, salvage3)
    @test round(c3; digits=1) == CostFactors(
        10.5, 1.1, 2.2, 3.3, 4.4, -0.6
    )
    @test round(c3; sigdigits=2) == CostFactors(
        11.0, 1.1, 2.2, 3.3, 4.4, -0.56
    )

end


@testset "Economics: entire Microgrid" begin
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
    # LCOE
    @test round(costs0.lcoe; digits=3) == 0.246 # $/kWh, without discount
    @test round(costs5.lcoe; digits=3) == 0.297 # $/kWh, with 5% discount
    # Cost factors of the system
    sys_c0M = CostFactors(
        total=41.697,
        investment=11.07,
        replacement=6.75,
        om=8.244,
        fuel=16.772,
        salvage=-1.139
    )
    sys_c5M = CostFactors(
        total=28.353,
        investment=11.07, # same as with no discount
        replacement=3.516,
        om=4.648,
        fuel=9.455,
        salvage=-0.336
    )
    @test round(costs0.system/1e6; digits=3) == sys_c0M # M$, without discount
    @test round(costs5.system/1e6; digits=3) == sys_c5M # M$, with 5% discount
end