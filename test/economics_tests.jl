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

    ctot = c1 + c2
    @test ctot.total == c1.total + c2.total
    @test ctot.investment == c1.investment + c2.investment
    @test ctot.replacement == c1.replacement + c2.replacement
    @test ctot.om == c1.om + c2.om
    @test ctot.fuel == c1.fuel + c2.fuel
    @test ctot.salvage == c1.salvage + c2.salvage
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
    # TODO: test LCOE
end