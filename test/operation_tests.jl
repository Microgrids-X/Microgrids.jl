# Tests for Microgrid operation simulation, in particular operation statistics

@testset "Microgrid operation" begin
    # Relative tolerance for comparing statistics:
    rtol = 1e-10

    # Parameters, similar to typical_parameters, with a few differences:
    # dt, fuel_intercept, charge/discharge_rate and ratings

    # Project
    lifetime = 25 # yr
    discount_rate = 0.05
    dt_1h = 1. # h
    dt_30m = 0.5 # h

    # Parameters common to all Components
    replacement_price_ratio = 1.0
    salvage_price_ratio = 1.0

    # Dispatchable generator
    power_rated_gen = 500.
    fuel_intercept_0 = 0.00 # constant efficiency generator
    fuel_intercept_5 = 0.05
    fuel_slope = 0.240
    fuel_price = 1.
    investment_price_gen = 400.
    om_price_gen = 0.02
    lifetime_gen = 15000.
    load_ratio_min = 0.0
    fuel_unit = "L"

    # Battery energy storage
    energy_rated = 6000.
    investment_price_sto = 350.
    om_price_sto = 10.
    lifetime_sto = 15.
    lifetime_cycles = 3000.
    charge_rate = 1.0
    charge_rate_slow = 1.0/30 # 200 kW charging
    discharge_rate = 1.0
    discharge_rate_slow = 1.0/30 # 200 kW discharging
    loss_factor_sto = 0.05
    SoC_min = 1/6 # 1000 kWh
    SoC_ini = 1/3 # 2000 kWh

    # Photovoltaic generation
    power_rated_pv = 1000. # set in x
    investment_price_pv = 1200.
    om_price_pv = 20.
    lifetime_pv = 25.
    derating_factor_pv = 1.

    ### Synthetic time series (see notebook for their design)
    Pload_max = 1000.0 # kW
    ndays = 7

    # slow time step
    K_1h = Int(ndays*24/dt_1h)

    t_1h = (0:K_1h-1)*dt_1h
    @assert length(t_1h) == K_1h
    td_1h = t_1h/24

    irradiance_1h = @. 0.5* (1 - cos(2*pi*t_1h/24)) # in [0,1]
    Pload_1h = collect(1 .- td_1h/ndays)*Pload_max

    # faster time step
    K_30m = Int(ndays*24/dt_30m)
    t_30m = (0:K_30m-1)*dt_30m
    @assert length(t_30m) == K_30m
    td_30m = t_30m/24

    irradiance_30m = @. 0.5* (1 - cos(2*pi*t_30m/24)) # in [0,1]
    Pload_30m = collect(1 .- td_30m/ndays)*Pload_max
    @assert length(irradiance_30m) == K_30m
    @assert length(Pload_30m) == K_30m

    ### Components:
    # generator with constant efficiency
    gen0 = DispatchableGenerator(power_rated_gen,
        fuel_intercept_0, fuel_slope, fuel_price,
        investment_price_gen, om_price_gen, lifetime_gen,
        load_ratio_min,
        replacement_price_ratio, salvage_price_ratio, fuel_unit)
    gen = DispatchableGenerator(power_rated_gen,
        fuel_intercept_5, fuel_slope, fuel_price,
        investment_price_gen, om_price_gen, lifetime_gen,
        load_ratio_min,
        replacement_price_ratio, salvage_price_ratio, fuel_unit)

    batt = Battery(energy_rated,
        investment_price_sto, om_price_sto, lifetime_sto, lifetime_cycles,
        charge_rate, discharge_rate, loss_factor_sto, SoC_min, SoC_ini,
        replacement_price_ratio, salvage_price_ratio)
    # battery with slow charging capability
    batt_slow_charge = Battery(energy_rated,
        investment_price_sto, om_price_sto, lifetime_sto, lifetime_cycles,
        charge_rate_slow, discharge_rate, loss_factor_sto, SoC_min, SoC_ini,
        replacement_price_ratio, salvage_price_ratio)
    # battery with slow discharging capability
    batt_slow_discharge = Battery(energy_rated,
        investment_price_sto, om_price_sto, lifetime_sto, lifetime_cycles,
        charge_rate, discharge_rate_slow, loss_factor_sto, SoC_min, SoC_ini,
        replacement_price_ratio, salvage_price_ratio)

    pv = Photovoltaic(power_rated_pv, irradiance_1h,
        investment_price_pv, om_price_pv,
        lifetime_pv, derating_factor_pv,
        replacement_price_ratio, salvage_price_ratio)
    # halved PV plant
    pv_half = Photovoltaic(power_rated_pv/2, irradiance_1h,
        investment_price_pv, om_price_pv,
        lifetime_pv, derating_factor_pv,
        replacement_price_ratio, salvage_price_ratio)
    pv_30m = Photovoltaic(power_rated_pv, irradiance_30m,
        investment_price_pv, om_price_pv,
        lifetime_pv, derating_factor_pv,
        replacement_price_ratio, salvage_price_ratio)

    # Microgrid variants
    proj = Project(lifetime, discount_rate, dt_1h, "€")
    proj_30m = Project(lifetime, discount_rate, dt_30m, "€")

    mg_base = Microgrid(proj, Pload_1h, gen, batt, [pv])
    # Microgrid with two halved PV plants → should be identical to base
    mg_2pv = Microgrid(proj, Pload_1h, gen, batt, [pv_half, pv_half])
    # Microgrid with 30 min time step → should give similar stats
    mg_30m = Microgrid(proj_30m, Pload_30m, gen, batt, [pv_30m])
    # Microgrid with constant efficiency generator → less fuel but same operating hours
    mg_gen0 = Microgrid(proj, Pload_1h, gen0, batt, [pv])
    # Microgrids with battery with slow (dis)-charging capability
    mg_slow_charge = Microgrid(proj, Pload_1h, gen, batt_slow_charge, [pv])
    mg_slow_discharge = Microgrid(proj, Pload_1h, gen, batt_slow_discharge, [pv])

    "simulate operation statistics, with optional smoothing of aggregation"
    function simulate_stats(mg, smoothing=NoSmoothing)
        # Run the microgrid operation
        oper_traj = operation(mg)
        # Aggregate the operation variables
        oper_stats = aggregation(mg, oper_traj, smoothing)
        return oper_stats
    end
    function simulate_Esto(mg)
        # Run the microgrid operation
        oper_traj = operation(mg)
        return oper_traj.Ebatt
    end

    "test equality of two structs, by each property"
    function test_properties(a,b; title="property", rtol)
        @testset "$title $name" for name in propertynames(a)
            va = getproperty(a, name)
            vb = getproperty(b, name)
            @test va ≈ vb rtol=rtol
        end
    end

    stats_base = OperationStats(
        # Load statistics
        79511.7634897933, 4988.236510206714, 477.01053219215305, 24.0, 10.0, 0.05903238473617412,
        # Dispatchable generator statistics
        17702.87093198049, 48.0, 5448.689023675317,
        # Energy storage (e.g. battery) statistics
        1.9836794636446038, 14497.180620961004, 9306.972942774242, 1190.2076781867618,
        # Non-dispatchable (typ. renewables) sources statistics
        17000.89976400046, 928.5714285714286, 0.20239166385714832,
        84000.0, 66999.10023599953, 0.7773553225963482
    )

    stats_base_show = """OperationStats with fields:
    - served_energy: 79512.0 kWh
    - shed_energy: 4988.2 kWh
    - shed_max: 477.01 kW
    - shed_hours: 24.0 h
    - shed_duration_max: 10.0 h
    - shed_rate: 0.059032 in [0,1]
    - gen_energy: 17703.0 kWh
    - gen_hours: 48.0 h
    - gen_fuel: 5448.7 L
    """ *
    "- storage_cycles: 1.9837 \n" *
    """- storage_char_energy: 14497.0 kWh
    - storage_dis_energy: 9307.0 kWh
    - storage_loss_energy: 1190.2 kWh
    - spilled_energy: 17001.0 kWh
    - spilled_max: 928.57 kW
    - spilled_rate: 0.20239 in [0,1]
    - renew_potential_energy: 84000.0 kWh
    - renew_energy: 66999.0 kWh
    - renew_rate: 0.77736 in [0,1]
    """

    @testset "Sanity check stats_base" begin
        expected_solar_energy = 0.5 * power_rated_pv * derating_factor_pv * K_1h*dt_1h
        @test stats_base.renew_potential_energy ≈ expected_solar_energy rtol=1e-3
    end

    @testset "Microgrid base case" begin
        stats_base_actual = simulate_stats(mg_base)
        test_properties(stats_base_actual, stats_base; title="OperStats", rtol=rtol)
        @testset "show OperationStats multiline" begin
            # `sprint` pastes the multi-line display in a String
            @test sprint(io -> show(io, "text/plain", stats_base_actual)) == stats_base_show
        end
        @testset "Esto" begin
            @test simulate_Esto(mg_base) ≈ [2000.0, 1000.0, 1000.0, 1000.0,
            1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0,
            1046.0171484396883, 1113.8742912968312, 1171.2009635460433,
            1186.7296970103184, 1126.7107571332563, 1000.0, 1000.0, 1000.0,
            1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0,
            1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1047.482863920753,
            1176.1068354802662, 1357.8382696342403, 1561.409698205669,
            1754.450656169167, 1905.693675347728, 1987.1051106970529,
            1974.6051106970529, 1841.735109375877, 1579.2351093758773,
            1187.1051106970544, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0,
            1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1079.1666666666663,
            1262.3638163017044, 1526.7020735755032, 1844.147793443763,
            2183.4335077294772, 2512.1887514072614, 2799.1460563001083,
            3016.2717773637182, 3140.676539268481, 3156.1751095017025,
            3043.675109501702, 2801.545110822879, 2439.045110822879,
            1974.064050699942, 1431.9007137131118, 1000.0, 1000.0, 1000.0,
            1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1094.6652368998873,
            1309.5461892808398, 1628.4576246301635, 2028.5101676182471,
            2481.6701732007928, 2956.6701732007928, 3421.139702592862,
            3843.8112931999945, 4196.6512999778915, 4456.77034759694,
            4607.983203544447, 4641.911774973018, 4549.781776294195,
            4337.281776294197, 4022.30071617126, 3630.13737918443,
            3191.776320382669, 2741.776320382669, 2315.9152615809076,
            1948.7519245940769, 1671.2708644711392, 1508.77086447114,
            1479.1408657923162, 1580.9265800780292, 1811.306102692202,
            2161.9013407874404, 2616.5270618510513, 3152.2938905534215,
            3741.1681818502525, 4351.882467564538, 4952.066282670894,
            5510.4521589923115, 5999.0064514844935, 6000.0, 6000.0, 6000.0,
            6000.0, 5937.5, 5772.5189398770635, 5530.355602890233,
            5241.994544088472, 4941.994544088472, 4666.133485286711,
            4448.970148299881, 4321.489088176942, 4308.989088176941,
            4417.895279848481, 4655.395279848481, 5021.489088176941,
            5507.798611986463, 6000.0, 6000.0, 6000.0, 6000.0, 6000.0, 6000.0,
            6000.0, 6000.0, 6000.0, 6000.0, 6000.0, 6000.0, 5985.018939877064,
            5892.855602890233, 5754.494544088473, 5604.494544088473,
            5478.633485286712, 5411.470148299881, 5431.844427236271,
            5556.249189141032, 5800.869666526858, 6000.0, 6000.0, 6000.0,
            6000.0, 6000.0, 6000.0, 6000.0, 6000.0, 6000.0, 6000.0, 6000.0,
            6000.0, 6000.0, 6000.0, 6000.0, 6000.0, 6000.0, 6000.0] rtol=rtol
        end
    end

    @testset "Microgrid with two halved PV plants" begin
        test_properties(simulate_stats(mg_2pv), stats_base; title="OperStats", rtol=rtol)
    end

    @testset "Microgrid with 30 min time step" begin
        # with timestep 1h → 30 min, the largest difference is on gen_hours: 48.0 → 45.5 h
        stats_30m = OperationStats( # difference with stats_base is at most 6% (e.g. gen_hours)
            # Load statistics
            79366.0717368912, 4883.92826310879, 477.01053219215316, 24.5, 10.0, 0.057969474933042026,
            # Dispatchable generator statistics
            17556.50106581271, 45.5, 5351.060255795051,
            # Energy storage (e.g. battery) statistics
            1.983957030705077, 14498.929293441983, 9308.55507501894, 1190.3742184230432,
            # Non-dispatchable (typ. renewables) sources statistics
            17000.05511049846, 928.5714285714286, 0.202381608458315,
            84000.0, 66999.94488950154, 0.7787908525444628
        )
        test_properties(simulate_stats(mg_30m), stats_base; title="OperStats", rtol=0.06)
        test_properties(simulate_stats(mg_30m), stats_30m; title="OperStats", rtol=rtol)
    end

    @testset "Microgrid with constant efficiency generator" begin
    # → expecting less fuel but same operating hours (and all other stats)
        stats_gen0 = OperationStats( # same as stats_base, but with less fuel (-1200 L)
            # Load statistics
            79511.7634897933, 4988.236510206714, 477.01053219215305, 24.0, 10.0, 0.05903238473617412,
            # Dispatchable generator statistics
            17702.87093198049, 48.0, 5448.689023675317 - power_rated_gen*fuel_intercept_5*48.0,
            # Energy storage (e.g. battery) statistics
            1.9836794636446038, 14497.180620961004, 9306.972942774242, 1190.2076781867618,
            # Non-dispatchable (typ. renewables) sources statistics
            17000.89976400046, 928.5714285714286, 0.20239166385714832,
            84000.0, 66999.10023599953, 0.7773553225963482
        )
        test_properties(simulate_stats(mg_gen0), stats_gen0; title="OperStats", rtol=rtol)
    end

    @testset "Microgrids with battery with slow charging capability" begin
        # small increase of load shedding and spilled energy
        stats_slow_charge = OperationStats(
            # Load statistics
            79451.41962426782, 5048.580375732201, 477.01053219215305, 25.0, 10.0, 0.059746513322274555,
            # Dispatchable generator statistics
            19360.131440266126, 54.0, 5996.431545663869,
            # Energy storage (e.g. battery) statistics
            1.6823453629759861, 12598.775786748713, 7589.368568963121, 1009.4072177855915,
            # Non-dispatchable (typ. renewables) sources statistics
            18899.304598212748, 799.9819620218456, 0.22499172140729462,
            84000.0, 65100.69540178725, 0.7563274321362443
        )
        test_properties(simulate_stats(mg_slow_charge), stats_slow_charge; title="OperStats", rtol=rtol)
    end

    @testset "Microgrids with battery with slow discharging capability" begin
        # As a side note, the *shedding rate is lower than the base case* (this relates to the undersized generator...)
        # → idea for an improved rule-based energy management...
        stats_slow_discharge = OperationStats(
            # Load statistics
            80171.36085007957, 4328.639149920452, 357.1428571428571, 22.0, 10.0, 0.05122649881562664,
            # Dispatchable generator statistics
            19764.429914552697, 67.0, 6418.463179492646,
            # Energy storage (e.g. battery) statistics
            1.737721284296193, 12947.644091066019, 7905.0113204883, 1042.6327705777185,
            # Non-dispatchable (typ. renewables) sources statistics
            18550.436293895444, 928.5714285714286, 0.2208385273082791,
            84000.0, 65449.563706104556, 0.7534726901853122
        )
        test_properties(simulate_stats(mg_slow_discharge), stats_slow_discharge; title="OperStats", rtol=rtol)
    end

    # Smoothed aggregation
    @testset "Smoothing of the discontinuous aggregated statistics" begin
        stats_relax_010 = OperationStats(
            # Load statistics: hours 24.0 → 23.28 h, duration max 10.0 → 9.32 h
            79511.7634897933, 4988.236510206714, 477.01053219215305, 23.281128854031188, 9.32435500456943, 0.05903238473617412,
            # Dispatchable generator statistics: hours 48 → 46.36 h, fuel 5448.68 → 5407.73 L
            17702.87093198049, 46.3618268038303, 5407.734693771075,
            # Energy storage (e.g. battery) statistics
            1.9836794636446038, 14497.180620961004, 9306.972942774242, 1190.2076781867618,
            # Non-dispatchable (typ. renewables) sources statistics
            17000.89976400046, 928.5714285714286, 0.20239166385714832,
            84000.0, 66999.10023599953, 0.7773553225963482
        )

        stats_relax_050 = OperationStats(
            # Load statistics: hours 24.0 → 16.496 h, duration max 10.0 → 7.59 h
            79511.7634897933, 4988.236510206714, 477.01053219215305, 16.49645098500673, 7.592538201845011, 0.05903238473617412,
            # Dispatchable generator statistics: hours 48 → 40.83 h, fuel 5448.68 → 5269.4 L
            17702.87093198049, 40.83038168444497, 5269.448565786442,
            # Energy storage (e.g. battery) statistics
            1.9836794636446038, 14497.180620961004, 9306.972942774242, 1190.2076781867618,
            # Non-dispatchable (typ. renewables) sources statistics
            17000.89976400046, 928.5714285714286, 0.20239166385714832,
            84000.0, 66999.10023599953, 0.7773553225963482
        )

        # And then set smoothing gain to 2.0:
        stats_relax_050_g2 = OperationStats(
            # Load statistics: hours 16.496 → 32.99 h, duration max 7.59 → 15.2 h (×2)
            79511.7634897933, 4988.236510206714, 477.01053219215305, 32.99290197001346, 15.185076403690022, 0.05903238473617412,
            # Dispatchable generator statistics: hours 40.83 → 81.66 h (×2), fuel 5269.4 → 6290.2 L
            17702.87093198049, 81.66076336888995, 6290.208107897567,
            # Energy storage (e.g. battery) statistics
            1.9836794636446038, 14497.180620961004, 9306.972942774242, 1190.2076781867618,
            # Non-dispatchable (typ. renewables) sources statistics
            17000.89976400046, 928.5714285714286, 0.20239166385714832,
            84000.0, 66999.10023599953, 0.7773553225963482
        )

        # on that simulation, relaxing with ε=1% has no effect
        # (explanation: relaxed variables are always above 1% of their range)
        test_properties(simulate_stats(mg_base, Smoothing(transition=0.01)),
                        stats_base; title="OperStats", rtol=rtol)
        test_properties(simulate_stats(mg_base, Smoothing(transition=0.10)),
                        stats_relax_010; title="OperStats", rtol=rtol)
        test_properties(simulate_stats(mg_base, Smoothing(transition=0.50)),
                        stats_relax_050; title="OperStats", rtol=rtol)
        test_properties(simulate_stats(mg_base, Smoothing(transition=0.50, gain=2.0)),
                        stats_relax_050_g2; title="OperStats", rtol=rtol)
    end
end