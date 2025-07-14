# Tests for Microgrid system description and project information structure

@testset "Microgrid project" begin

    proj = Project(25, 0.05, 1.0, "€") # LinearSalvage default
    @test proj.salvage_type == LinearSalvage

    # Test mutability and copy of Porject
    proj1 = copy(proj)
    proj1.discount_rate = 0.20
    proj1.currency = "EUR"
    @test proj.discount_rate == 0.05
    @test proj.currency == "€"

end

@testset "Microgrid system description" begin
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
    proj = Project(lifetime, 0.0, timestep, "€") # LinearSalvage default

    mg = Microgrid(proj, Pload, generator, battery, [photovoltaic])

    mg_shallow = copy(mg)
    mg_shallow.project.discount_rate = 0.20
    # Original Project and thus Microgrid is changed:
    @test proj.discount_rate == 0.20
    @test mg.project.discount_rate == 0.20

    mg_shallow.storage.energy_rated = 2*energy_rated
    # Original Battery and thus Microgrid.storage is changed:
    @test battery.energy_rated == 2*energy_rated
    @test mg.storage.energy_rated == 2*energy_rated
    # Now copy Battery, so that the original one is not changed:
    mg_shallow.storage = copy(mg_shallow.storage)
    mg_shallow.storage.energy_rated = 3*energy_rated
    @test mg_shallow.storage.energy_rated == 3*energy_rated
    @test mg.storage.energy_rated == 2*energy_rated # unchanged

    # Same for non-dispatchable sources:
    mg_shallow.nondispatchables[1].power_rated = 2*power_rated_pv
    # Original Photovoltaic and thus Microgrid.nondispatchables is changed:
    @test photovoltaic.power_rated == 2*power_rated_pv
    @test mg.nondispatchables[1].power_rated == 2*power_rated_pv
    # Now copy both the Vector and the PV
    mg_shallow.nondispatchables = [copy(photovoltaic)]
    mg_shallow.nondispatchables[1].power_rated = 3*power_rated_pv
    @test mg_shallow.nondispatchables[1].power_rated == 3*power_rated_pv
    @test photovoltaic.power_rated == 2*power_rated_pv # unchanged

end