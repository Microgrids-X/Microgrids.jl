# Tests of Microgrids.jl with respect to optimization (e.g. differientation)
using ForwardDiff

@testset "Differentiate simulation(mg)" begin
    include("typical_parameters.jl")

    # Dummy time series (only 3 hours of operation, so that cost value below is meaningless)
    Pload = [1., 1., 1.]
    IT = [0, 0.5, 1.0]

    """simulate microgrid of size `x` (size=3) and returns its Net Present Cost"""
    function sim_npc(x)
        power_rated_DG = x[1]
        energy_max = x[2]
        power_rated_PV = x[3]

        # Components:
        dieselgenerator = DieselGenerator(power_rated_DG, min_load_ratio, F0, F1, fuel_cost, investiment_cost_DG, om_cost_DG, replacement_cost_DG, salvage_cost_DG, lifetime_DG)
        battery = Battery(energy_initial, energy_max, energy_min, power_min, power_max, loss, investiment_cost_BT, om_cost_BT, replacement_cost_BT, salvage_cost_BT, lifetime_BT, lifetime_thrpt)
        photovoltaic = Photovoltaic(power_rated_PV, fPV, IT, IS, investiment_cost_PV, om_cost_PV, replacement_cost_PV, salvage_cost_PV, lifetime_PV)

        # Microgrid project, with discount rate:
        proj = Project(lifetime, discount_rate, timestep)
        mg = Microgrid(proj, Pload, dieselgenerator, battery, [photovoltaic])

        results = simulate(mg)
        return results.costs.npc
    end

    x = [power_rated_DG, energy_max, power_rated_PV]

    # Just a appetizer to check that cost computation works
    @test round(sim_npc(x)/1e6; digits=3) == 15.023 # M$

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
    @test round.(grad_approx; digits=3) == [282.358, 624.843, 1481.879] # $/kW or $/kWh

    # Now gradient computation with ForwardDiff:
    grad_fd = ForwardDiff.gradient(sim_npc, x)
    @test round.(grad_fd; digits=3) == [282.358, 624.843, 1481.879] # $/kW or $/kWh

end