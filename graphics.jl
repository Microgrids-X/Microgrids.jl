 td = collect((0:nsteps-1)/24);

"""
This file must be included after the instructions:
 - Using(Microgrids.jl)
 - Include ("../examples/data/Microgrid_Wind-Solar-H2_data.jl")
"""

"""
Plot the instantaneous power sharing of a microgrid furing one year
plot_oper_traj(mg, oper_traj)
"""
function plot_oper_traj(mg, oper_traj)
    
   
    
    fig,( ax1, ax2) = plt.subplots(2,1,sharex=true)

    y1 = min.(oper_traj.Prenew_pot,mg.load)
    y2 = max.( oper_traj.Pbatt,0.0)
    y3 = oper_traj.Pfc
    y4 = abs.(min.( oper_traj.Pbatt,0.0))
    y5 = oper_traj.Pelyz
    y6 = oper_traj.power_shedding
    y7 = oper_traj.power_curtailment
    
 
    y=np.vstack([y1,y2,y3,y4,y5,y6,y7])
    fieldNames =(["Pren used by the load","battery_discharge","fuel cell","battery_charge", "elyz","shedding","spillage"])
    
    fieldColors = (["orange","mediumseagreen","hotpink","seagreen","orchid","black","red"])
         
    
    ax1.plot(td, oper_traj.Prenew_pot,lw=2, "tab:orange", label="ren")
    ax1.plot(td, mg.load,lw=2, label="load")
    stacks = ax1.stackplot(td,y,labels=fieldNames, colors = fieldColors,alpha=0.8)
   
    stacks[2].set_hatch("--")
    stacks[4].set_hatch("++")
    stacks[6].set_hatch("\\")   
    stacks[7].set_hatch("\\") 
    
    ax3=ax1.twinx()
    ax3.set_ylim(0, 2)
    ax3.plot(td, oper_traj.Ebatt[1:8760]/mg.storage.energy_rated, ls="--","C2",label="hourly battery level")
    
 #  ax4=ax2.twinx()
   
    ax2.plot(td, oper_traj.LoH[1:8760], "C2",label="level of H2")
    
  #  ax2.stairs(z,tm,fill="True",label="level of H2 at the end of each month")
         
    
    ax1.set_title("hourly power repartition")
    ax2.set_title("H2 level")
    ax1.legend()
    ax1.grid(true)
    ax1.set(ylabel="kW")
    ax2.legend()
    ax3.legend()
    ax2.grid(true)
    ax2.set(ylabel="kg")
    


       fig.tight_layout()  
       pygui(true)
       plt.show()
     return fig, ax1, ax2
    
end

function comp_oper_traj(mg1,mg2, oper_traj1, oper_traj2)
 
    ns=Int64(nsteps/730)
       z=zeros(Float64,ns)
  tm=zeros(Float64,ns +1)
      tm[end]=nsteps/24
   for i=1:ns
      z[i]= oper_traj1.LoH[i*730+1]
      tm[i]=(i-1)*730/24 
  end
      fig,( ax1, ax2) = plt.subplots(2,1,sharex=true)
      y1 = min.(oper_traj1.Prenew_pot,mg1.load)
      y2 = max.( oper_traj1.Pbatt,0.0)
      y3 = oper_traj1.Pfc
      y4 = abs.(min.( oper_traj1.Pbatt,0.0))
      y5 = oper_traj1.Pelyz
      y6 = oper_traj1.power_shedding
      y7 = oper_traj1.power_curtailment
      
      Y1=np.vstack([y1,y2,y3,y4,y5,y6,y7])
      
   
      y1 = min.(oper_traj2.Prenew_pot,mg2.load)
      y2 = max.( oper_traj2.Pbatt,0.0)
      y3 = oper_traj2.Pfc
      y4 = abs.(min.( oper_traj2.Pbatt,0.0))
      y5 = oper_traj2.Pelyz
      y6 = oper_traj2.power_shedding
      y7 = oper_traj2.power_curtailment
      
      Y2=np.vstack([y1,y2,y3,y4,y5,y6,y7])
   
          
          fieldNames =(["Pren used by the load","battery_discharge","fuel cell","battery_charge", "elyz","shedding","spillage"])
      
      fieldColors = (["orange","green","pink","green","magenta","black","red"])
           
      
      ax1.plot(td, oper_traj1.Prenew_pot,lw=2, "tab:orange", label="ren")
      ax1.plot(td, mg1.load,lw=2, label="load")
      stacks1 = ax1.stackplot(td,Y1,labels=fieldNames, colors = fieldColors,alpha=0.8)
       ax2.plot(td, oper_traj2.Prenew_pot,lw=2, "tab:orange", label="ren")
      ax2.plot(td, mg2.load,lw=2, label="load")
      stacks2 = ax2.stackplot(td,Y2,labels=fieldNames, colors = fieldColors,alpha=0.8)
      
      
      stacks1[2].set_hatch("-")
      stacks1[4].set_hatch("+")
      stacks1[6].set_hatch("\\")   
      stacks1[7].set_hatch("\\") 
      
      stacks2[2].set_hatch("-")
      stacks2[4].set_hatch("+")
      stacks2[6].set_hatch("\\")   
      stacks2[7].set_hatch("\\") 
      
      ax3=ax1.twinx()
      ax3.set_ylim(0, 2)
      ax3.plot(td, oper_traj1.Ebatt[1:end-1]/mg1.storage.energy_rated, ls="--","C2",label="hourly battery level")
  
      
      
      ax4=ax2.twinx()
      ax4.set_ylim(0, 2)
      ax4.plot(td, oper_traj2.Ebatt[1:end-1]/mg2.storage.energy_rated, ls="--","C2",label="hourly battery level")
  
      
             
      
      ax1.set_title("hourly power repartition")
      ax2.set_title("hourly power repartition")
      ax1.legend()
      ax1.grid(true)
      ax1.set(ylabel="kW")
      ax2.legend()
      ax3.legend()
       ax4.legend()
      ax2.grid(true)
      ax2.set(ylabel="kw")
       fig.tight_layout()
     
       return fig,(ax1, ax2,ax3,ax4)
      
  end
  """
    plot_cashflow(cashflow)

TBW
"""
function plot_cashflow(cashflow)
    fig,ax=plt.subplots()
    
    years=collect(range(0,25,26))
    x=[cashflow.pv, cashflow.wind, cashflow.fc ,cashflow.elyz, cashflow.batt, cashflow.h2_tank]
    names=["PV" "Wind" "fc" "elyz" "batt" "h2_tank"]
    
    
        bottom = zeros(Float64,26)
    
        
        for i=1:6
        p = ax.bar(years,x[i], 0.8, label=names[i], bottom=bottom,align="edge",capstyle="projecting")
        bottom .+= x[i]
        end
        ax.set_title("cashflow diagrams")
    ax.legend(loc="lower right")
      ax.set_ylim(-3e7,2e7)
        pygui(true)
       plt.show()
        return fig, ax
    end

    function plot_cashflow_comp(cashflow_1,cashflow_2)
        fig,ax=plt.subplots()
        
        years=collect(range(0,25,26))
        x1=[cashflow_1.wind, cashflow_1.pv, cashflow_1.fc ,cashflow_1.elyz, cashflow_1.batt, cashflow_1.h2_tank]
        x2=[cashflow_2.wind, cashflow_2.pv, cashflow_2.fc ,cashflow_2.elyz, cashflow_2.batt, cashflow_2.h2_tank]
        names=["Wind" "PV" "fc" "elyz" "batt" "h2_tank"]
        
        width = 0.4  # the width of the bars
        color=["lavender","orange","darkviolet","lightcyan","darkgreen","cornflowerblue"]
            bottom_1 = zeros(Float64,26)
             bottom_2 = zeros(Float64,26)
            
            for i=1:6
            rects_1 = ax.bar(years,x1[i],width, label=names[i], bottom=bottom_1,color=color[i],align="edge",capstyle="projecting",ec="darkgreen")
            rects_2 = ax.bar(years .+ width,x2[i],width, label=names[i], bottom=bottom_2,color=color[i],align="edge",capstyle="projecting",ec="darkred")
            bottom_1 .+= x1[i]
            bottom_2 .+= x2[i]
            #ax.bar_label(rects, padding=3)
            end
            ax.set_title("cashflow diagrams")
        ax.legend(loc="lower right")
            pygui(true)
           plt.show()
            plt.savefig("cash_flow_comp.png")
            return fig, ax
        end