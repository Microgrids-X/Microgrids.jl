using Microgrids
using NLopt


"Multi-objective criterion for microgrid performance: lcoe, shedding rate"
function obj_multi(x= X,capex=capex_def,dispatch=dispatch_1,load=Pload,cf_pv=Ppv1k,cf_wind=cf_wind)
         mg=new_microgrid(Sizing((1000*x)...),capex,ini_filling_state,load,cf_wind,cf_pv)
    mg,traj,stats, costs = simulate_double(mg,capex,dispatch)
    # Extract KPIs of interest
    lcoe = costs.lcoe # $/kWh
    shed_rate = stats.shed_rate; # in [0,1]
    return lcoe, shed_rate
end

"""Mono-objective criterion: LCOE + penalty if shedding rate > `shed_max`

with signature adapted to NLopt with `grad` as 2nd argument

load shedding penalty threshold `shed_max` should be in  [0,1[
"""
function obj(x,grad, shed_max,capex=capex_def,dispatch=dispatch_1, w_shed_max=1e5,load=Pload,cf_pv=Ppv1k,cf_wind=cf_wind)
    lcoe, shed_rate =  obj_multi(x,capex,dispatch,load,cf_pv,cf_wind)
    over_shed = shed_rate - shed_max
    if over_shed > 0.0
        penalty = w_shed_max*over_shed
    else
        penalty = 0.0
    end
    J = lcoe + penalty
end

"""Optimize sizing of microgrid based on the `obj` function

Parameters:
- `x0`: initial sizing (for the algorithms which need them)
- `shed_max`: load shedding penalty threshold (same as in `obj`)
- `algo` could be one of LN_SBPLX, GN_DIRECT, GN_ESCH...
- `maxeval`: maximum allowed number of calls to the objective function,
  that is to the microgrid simulation
- `xtol_rel`: termination condition based on relative change of sizing, see NLopt doc.
- `srand`: random number generation seed (for algorithms which use some stochastic search)

"""
function optim_mg(x0, shed_max, xmin,xmax, algo=:GN_CRS2_LM, maxeval=1000, xtol_rel=1e-8, srand=1,load=Pload,cf_pv=Ppv1k,cf_wind=cf_wind)
    nx = length(x0) # number of optim variables
    opt = Opt(algo, nx)
    NLopt.srand(srand)
    
    opt.lower_bounds = xmin
    opt.upper_bounds = xmax
    opt.min_objective = (x, grad) -> obj(x, grad, shed_max,capex_def,Microgrids.dispatch_1,1e5,load,cf_pv,cf_wind)
    opt.xtol_rel = xtol_rel
    opt.maxeval = maxeval
    
    (fopt, xopt, ret) = optimize(opt, x0)
    return xopt, ret, opt.numevals
end