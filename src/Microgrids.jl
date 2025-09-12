module Microgrids
    using CSV, DataFrames
    using DocStringExtensions



    @template (TYPES) =
    """
    $(TYPEDEF)
    $(DOCSTRING)

    ### Parameters: 
    $(FIELDS)
    """

    include("components.jl")
    include("dispatch.jl")
    include("production.jl")
    include("operation.jl")
    include("economics.jl")
    include("../examples/data/Microgrid_Wind-Solar-H2_data.jl")


    import Base.@kwdef # backport Julia 1.9 syntax to 1.6-1.8 versions

    export simulate,save_mg,load_mg,new_microgrid,simulate_pnl,
        NonDispatchableSource, ProductionUnit, Tank, TankCompound, 
        Project, Sizing,DispatchableCompound,Battery, Photovoltaic, PVInverter, WindPower, Microgrid,
        capacity_from_wind,
        OperationTraj, OperationStats, increment,
        operation, aggregation, dispatch_1,dispatch_2, production, 
        CostFactors, MicrogridCosts, component_costs, economics, MicrogridCashflow
        

"""
    Create a Microgrid of size `x`
    with x=[energy_rated_sto, power_rated_pv, power_rated_wind, power_rated_elyz, power_rated_fc, capacity_rated_hy_tank]
    You can also specify capex prices and initial filling rate of either Hydrogen tank or batteries.
    new_microgrid(x::Vector{Float64}= X,capex::Vector{Float64}=capex_def,initial_fill_rate::Vector{Float64})

"""
        
function new_microgrid(sizing::Sizing = default_sizing,capex::Vector{Float64}=capex_def,initial_fill_rate::Vector{Float64}=ini_filling_state,load::Vector{Float64}=Pload,
    wind_speed::Vector{Float64}=wind_speed,irradiance::Vector{Float64}=irradiance)
    cf_wind = capacity_from_wind.(wind_speed; TSP=TSP_D52, Cp=Cp_D52, v_out=v_out, α=α_D52)
    project = Project(lifetime, discount_rate, timestep, "€")
  gen = ProductionUnit(sizing.Cgen,
      fuel_intercept, fuel_slope, fuel_price,
      capex[1], om_price_hour_gen, om_price_gen, lifetime_gen_y,lifetime_gen_h,lifetime_gen_starts,
      load_ratio_min_gen, replacement_price_ratio,
      salvage_price_ratio, input_unit_gen,output_unit_gen)
  ftank = Tank(sizing.Ftank,capex[2], om_price_ftank,lifetime_ftank,loss_factor_ftank,initial_fill_rate[1],
          fuel_min_ratio, fuel_max_ratio,fuel_price, replacement_price_ratio, salvage_price_ratio)
  fuel_cell = ProductionUnit(sizing.Cfc,cons_intercept_fc, cons_rate_fc,cons_price_fc,capex[3], om_price_hour_fc, om_price_fc,lifetime_fc_y,lifetime_fc_h,lifetime_fc_starts,
              load_min_ratio_fc,replacement_price_ratio, salvage_price_ratio,input_unit_fc,output_unit_fc)
  hytank = Tank(sizing.Htank,capex[4], om_price_hytank,lifetime_hytank,loss_factor_hytank,initial_fill_rate[2],
              LoH_min_ratio, LoH_max_ratio,hy_price,replacement_price_ratio, salvage_price_ratio)
  dispatchables = DispatchableCompound{Float64}([gen], [fuel_cell])
     tanks = TankCompound{Float64}(ftank,hytank)

  elyz = ProductionUnit(sizing.Cel,cons_intercept_elyz,cons_rate_elyz,cons_price_elyz, capex[5], om_price_hour_elyz,om_price_elyz, lifetime_elyz_y,lifetime_elyz_h,lifetime_elyz_starts,
  load_min_ratio_elyz,replacement_price_ratio, salvage_price_ratio,input_unit_elyz,output_unit_elyz)

  batt = Battery(sizing.Cbatt,
      capex[6], om_price_sto, lifetime_sto, lifetime_cycles,
      charge_rate, discharge_rate, loss_factor_sto, SoC_min, SoC_max,initial_fill_rate[3],
      replacement_price_ratio, salvage_price_ratio)
  pv = Photovoltaic(sizing.Cpv, irradiance,
      capex[7], om_price_pv,
      lifetime_pv, derating_factor_pv,
      replacement_price_ratio, salvage_price_ratio)
  windgen = WindPower(sizing.Cwind, cf_wind,
      capex[8], om_price_wind,
      lifetime_wind,
      replacement_price_ratio, salvage_price_ratio)

      hb=ProductionUnit(sizing.Hb,cons_intercept_hb, cons_rate_hb,cons_price_hb,capex[3], om_price_hour_hb, om_price_hb,lifetime_hb_y,lifetime_hb_h,lifetime_hb_starts,
      load_min_ratio_hb,replacement_price_ratio, salvage_price_ratio,input_unit_hb,output_unit_hb)


    mg = Microgrid(project, load,dispatchables,
[elyz,],hb,tanks,batt, [
    pv,
    windgen
    ])
   

  return mg
  end



    """
        simulate(mg::Microgrid, ε::Real=0.0)

    Simulate the technical and economic performance of a Microgrid `mg`.

    Discontinuous computations can optionally be relaxed (smoothed)
    using the relaxation parameter `ε`:
    - 0.0 means no relaxation (default value)
    - 1.0 yields the strongest relaxation

    Using relaxation (`ε` > 0) is recommended when using gradient-based optimization
    and then a “small enough” value between 0.05 and 0.30 is suggested.

    Returns:
    - Operational trajectories from `sim_operation` (should be optional in future version)
    - Operational statistics from `sim_operation`
    - Microgrid project costs from `sim_economics`
    """
    function simulate(mg::Microgrid, dispatch::Function,ε::Real=0.0)
        # Run the microgrid operation
        oper_traj = operation(mg, dispatch)
        # Aggregate the operation variables
        oper_stats = aggregation(mg, oper_traj, ε)

        # Eval the microgrid costs
        mg_costs = economics(mg, oper_stats)

        return oper_traj, oper_stats, mg_costs
    end
    function simulate_pnl(mg::Microgrid, dispatch::Function,Pnl_request::Vector{Float64},ε::Real=0.0)
        # Run the microgrid operation
        oper_traj = operation_pnl(mg, dispatch,Pnl_request)
        # Aggregate the operation variables
        oper_stats = aggregation(mg, oper_traj, ε)

        # Eval the microgrid costs
        mg_costs = economics(mg, oper_stats)

        return oper_traj, oper_stats, mg_costs
    end
    
# to update
    function save_mg(mg::Microgrid , stats::OperationStats, traj::OperationTraj, path::String,mg_name::String)
        mkdir(mg_name)
        proj_path=path*mg_name*"/"
        nsteps=Int(length(mg.load)/mg.project.timestep)
        size_frame=DataFrame(Cgen=Float64[],
        Cbatt=Float64[],Cpv=Float64[],Cwind=Float64[],Cfc=Float64[],Cel=Float64[],Htank=Float64[],Ftank=Float64[],Hb=Float64[])
        push!(size_frame,(mg.dispatchables.generator[1].power_rated,mg.storage.energy_rated,mg.nondispatchables[1].power_rated,mg.nondispatchables[2].power_rated,
        mg.dispatchables.fuel_cell[1].power_rated,mg.electrolyzer[1].power_rated,mg.tanks.h2Tank.capacity,mg.tanks.fuelTank.capacity,0.0))
        CSV.write(proj_path*"sizing.csv",size_frame)

        op_frame=DataFrame([name => [] for name in propertynames(stats)])
        push!(op_frame,collect(getproperty(stats, Symbol(name)) for name in propertynames(stats) ))
        CSV.write(proj_path*"stats.csv",op_frame)
        
        traj_frame=DataFrame([name => getproperty(traj,name)[1:nsteps] for name in propertynames(traj)])
        CSV.write(proj_path*"traj.csv",traj_frame)
        
    end

    function load_mg(project_path :: String)
        sizing_df = DataFrame(CSV.File(project_path*"/sizing.csv"))
        sizing= Sizing(sizing_df[1,1],sizing_df[1,2],sizing_df[1,3],sizing_df[1,4],sizing_df[1,5],sizing_df[1,6],sizing_df[1,7],sizing_df[1,8],0.0)
        
        traj_df = DataFrame(CSV.File(project_path*"/traj.csv"))
        traj= OperationTraj(traj_df[:,1],traj_df[:,2],traj_df[:,3],traj_df[:,4],traj_df[:,5],traj_df[:,6],zeros(8760),traj_df[:,7],traj_df[:,8],traj_df[:,9],traj_df[:,10],traj_df[:,11],traj_df[:,12],)
       
        if sizing.Ftank>0.0
            a = traj.LoF[1]/sizing.Ftank
        else
            a=0.0
        end
       
        if sizing.Htank>0.0
             b = traj.LoH[1]/sizing.Htank
        else
            b=0.0
        end
        if sizing.Cbatt>0.0
            c = traj.Ebatt[1]/sizing.Cbatt
       else
           c=0.0
       end
        mg=new_microgrid(sizing,capex_def,[a,b,c])
        stats = aggregation(mg, traj, 0.0)

        return mg,traj,stats

    end
end