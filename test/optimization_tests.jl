# Tests of Microgrids.jl with respect to optimization (e.g. differientation)
# with optional relaxation of discontinuous statistics
using ForwardDiff

@testset "Differentiate simulation(mg)" begin
    include("typical_parameters.jl")

    # Dummy time series (only 3 hours of operation, so that cost value below is meaningless)
    Pload = [1., 1., 1.]
    irradiance = [0, 0.5, 1.0]

    """simulate microgrid of size `x` (size=3) and returns its Net Present Cost

    Optionally smooth discontinuous statistics with relaxation parameter `ε`
    """
    function sim_npc(x, ε=0.0)
        power_rated_gen = x[1]
        energy_rated = x[2]
        power_rated_pv = x[3]

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

        # Microgrid project, with discount rate:
        proj = Project(lifetime, discount_rate, timestep, "€")
        mg = Microgrid(proj, Pload, generator, battery, [photovoltaic])

        oper_traj, oper_stats, mg_costs = simulate(mg, ε)
        return mg_costs.npc
    end

    x = [power_rated_gen, energy_rated, power_rated_pv]

    # Just a appetizer to check that cost computation works
    npc_expected = 15.023 # M$
    npc_expected_010 = 15.022 # M$
    @test round(sim_npc(x)/1e6; digits=3) == npc_expected
    @test round(sim_npc(x, 0.10)/1e6; digits=3) == npc_expected_010

    # Centered finite difference approximation of the gradient:
    dx = x*1e-3
    grad_approx = zero(x)
    for i in eachindex(x)
        x_right = copy(x)
        x_right[i] += dx[i]
        x_left = copy(x)
        x_left[i] -= dx[i]
        grad_approx[i] = (sim_npc(x_right) - sim_npc(x_left))/(2*dx[i])
    end

    # Test the finite difference approx of the gradient against recorded value
    grad_expected = [282.358, 624.843, 1481.879]
    grad_expected_010 = [281.879, 624.843, 1481.879] # notice that the gradient change is only with respect to power_rated_gen
    @test round.(grad_approx; digits=3) == grad_expected # $/kW or $/kWh

    # Now gradient computation with ForwardDiff:
    grad_ad = ForwardDiff.gradient(sim_npc, x)
    @test round.(grad_ad; digits=3) == grad_expected # $/kW or $/kWh
    grad_ad_010 = ForwardDiff.gradient(x -> sim_npc(x, 0.10), x)
    @test round.(grad_ad_010; digits=3) == grad_expected_010 # $/kW or $/kWh

end