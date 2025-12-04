# Load Microgrid prior using this file
using Interpolations, Base.Threads, LinearAlgebra, Base.Iterators, Distributions,FreqTables,StatsBase
using Microgrids

# Structure describing the state of the Microgrid
struct State
    Pnl::Float32
    LoH::Float32
    Ebatt::Float32
    δ::Int32
    # Custom constructor
    function State(x::Float32, y::Float32,s::Float32,z::Int32)
        return new(x, y,s,z)  # No allocations, GPU-safe
    end
end
"""
    State(a::NTuple{6, Float32})

TBW

function State(a::NTuple{4,Topt}){Topt<:Real}
    State(a[1], a[2], a[3], a[4])
end
"""

struct Time
    hi::Float32
    di::Int32
    mi::Int32
end

"""
    time_dyn_forward(time::Array{Float32, 1},dt::Float32,mi_l::Vector{Int32})


"""
function time_dyn_forward(time::Array{Float32,1}, dt::Float32, mi_l::Vector{Int32})

    if time[1] < 23.0
        next_hi = time[1] + dt
        next_di = time[2]
        next_mi = time[3]
    else
        next_hi = 0.0f0
        if time[2] < mi_l[Int(time[3])]
            next_di = time[2] + 1
            next_mi = time[3]
        else
            next_di = 1
            next_mi = time[3] + 1
        end
    end

    time[1] = next_hi
    time[2] = next_di
    time[3] = next_mi
end




"""
    time_dyn_backward!(time::Array{Float32, 1},dt::Float32,mi_l::Vector{Int32})

TBW
"""
function time_dyn_backward!(time::Array{Float32,1}, dt::Float32, mi_l::Vector{Int32})

    if time[1] > 0
        next_hi = time[1] - dt
        next_di = time[2]
        next_mi = time[3]
    else
        next_hi = 23.0
        if time[2] > 1
            next_di = time[2] - 1
            next_mi = time[3]
        else
            next_di = mi_l[Int(time[3])-1]
            next_mi = time[3] - 1
        end
    end

    time[1] = next_hi
    time[2] = next_di
    time[3] = next_mi
end



"""
    bt_bounds( x::State,mg::Microgrid)

TBW
"""
function bt_bounds(x, bt,dt)
   
        Pdis_max = min((x.Ebatt - bt.energy_rated * bt.SoC_min) / ((1 + bt.loss_factor) * dt), bt.energy_rated * bt.discharge_rate)
        Pcha_max = max((x.Ebatt - bt.energy_rated * bt.SoC_max) / ((1 - bt.loss_factor) * dt), -bt.energy_rated * bt.charge_rate)

    return Pcha_max, Pdis_max
end



"""
    h2_bounds( x::State,mg::Microgrid)

TBW
"""
function h2_bounds(x, fc,el,hytank,dt)

    
    Pel_min, Pfc_min = -el.power_rated * el.minimum_load_ratio, fc.power_rated * fc.minimum_load_ratio

   
    Pel_max = ifelse((x.LoH - hytank.capacity * hytank.max_filling_ratio) * el.consumption_slope / dt > Pel_min, 0.0f0, max((x.LoH - hytank.capacity * hytank.max_filling_ratio) * el.consumption_slope / dt, -el.power_rated))
 

   
    Pfc_max = ifelse((x.LoH - hytank.capacity * hytank.min_filling_ratio) / (fc.consumption_slope * dt) < Pfc_min, 0.0f0, min((x.LoH - hytank.capacity * hytank.min_filling_ratio) / (fc.consumption_slope * dt), fc.power_rated))
    

    return Pel_max, Pfc_max, Pel_min, Pfc_min
end

"""
    u_bounds( x::State,mg)

TBW
"""
function u_bounds(x, fc,el,hytank,bt,dt)
   
    Pel_max, Pfc_max, Pel_min, Pfc_min = h2_bounds(x, fc,el,hytank,dt)
    Pcha_max, Pdis_max = bt_bounds(x, bt,dt)

    if x.Pnl >= 0.0f0
        u_max = min(Pdis_max, x.Pnl)
        u_min = ifelse(Pfc_max >= x.Pnl, max(Pcha_max, x.Pnl - Pfc_max), min(Pdis_max, x.Pnl - Pfc_max))
    else
        u_max = max(min(x.Pnl - Pel_max, 0.0f0), Pcha_max)
        u_min = max(Pcha_max, x.Pnl - Pfc_max)
    end

    return u_min, u_max, Pel_min, Pfc_min
end


"""
    u_range!(u_range,x::State,mg,step=100,db=10,da=10)

TBW
"""
function u_range!(u_range, x, fc,el,hytank,bt, dt,step=100, db=10, da=10)
   
    u_min, u_max, Pel_min, Pfc_min = u_bounds(x, fc,el,hytank,bt,dt)

    na = da
    nb = db
    p1 = x.Pnl - Pfc_min
    p2 = x.Pnl - Pel_min
    nu = 0
    #  u_range=zeros(25)

    if x.Pnl >= 0.0f0
        if u_min != u_max

            if p1 <= u_max && p1 >= u_min
                if p1 >= 0.0f0

                    if 0.0f0 <= u_max && 0.0f0 >= u_min
                        na = min(Int(div(-u_min, step) + 2), da)
                        nb = min(Int(div(p1, step) + 2), db)
                        for i = 0:na-1
                            u_range[i+1] = u_min + (i * (-u_min)) / (na - 1)
                        end
                        for i = na:na+nb-2
                            u_range[i+1] = ((i - na + 1) * (p1)) / (nb - 1)
                        end
                        nu = na + nb - 1
                    else
                        nb = min(Int(div(p1 - u_min, step) + 2), db)
                        for i = 0:nb-1
                            u_range[i+1] = u_min + (i * (p1 - u_min)) / (nb - 1)
                        end
                        nu = nb
                    end
                else
                    na = min(Int(div(p1 - u_min, step) + 2), da)
                    for i = 0:na-1
                        u_range[i+1] = u_min + (i * (p1 - u_min)) / (na - 1)
                    end
                    nu = na
                end
            else
                if p1 > u_max
                    if 0.0f0 <= u_max && 0.0f0 >= u_min
                        na = min(Int(div(-u_min, step) + 2), da)
                        nb = min(Int(div(u_max, step) + 2), db)
                        for i = 0:na-1
                            u_range[i+1] = u_min + (i * (-u_min)) / (na - 1)
                        end
                        for i = na:na+nb-2
                            u_range[i+1] = ((i - na + 1) * (u_max)) / (nb - 1)
                        end
                        nu = na + nb - 1
                    else
                        nb = min(Int(div(u_max - u_min, step) + 2), db)
                        for i = 0:nb-1
                            u_range[i+1] = u_min + (i * (u_max - u_min)) / (nb - 1)
                        end
                        nu = nb

                    end
                end

            end
            if u_max == x.Pnl
                u_range[nu+1] = x.Pnl
                nu += 1
            end
        else
            u_range[1] = u_max
            nu = 1
        end

    else
        if u_min == u_max
            u_range[1] = u_max
            nu = 1
        else
            if p2 <= u_max && p2 >= u_min
                if x.Pnl >= u_min
                    if p1 >= u_min
                        na = min(Int(div(p1 - u_min, step) + 2), da)
                        nb = min(Int(div(u_max - p2, step) + 2), db)
                        for i = 0:na-1
                            u_range[i+1] = u_min + (i * (p1 - u_min)) / (na - 1)
                        end
                        u_range[na+1] = x.Pnl
                        for i = na+1:na+nb
                            u_range[i+1] = p2 + ((i - na - 1) * (u_max - p2)) / (nb - 1)
                        end
                        nu = na + nb + 1

                    else
                        nb = min(Int(div(u_max - p2, step) + 2), db)
                        u_range[1] = x.Pnl
                        for i = 0:nb-1
                            u_range[i+2] = p2 + (i * (u_max - p2)) / (nb - 1)
                        end
                        nu = nb + 1

                    end
                else
                    nb = min(Int(div(u_max - p2, step) + 2), db)
                    for i = 0:nb-1
                        u_range[i+1] = p2 + (i * (u_max - p2)) / (nb - 1)
                    end
                    nu = nb
                end
            else
                if x.Pnl <= u_max && x.Pnl >= u_min
                    if p1 >= u_min
                        na = min(Int(div(p1 - u_min, step) + 2), da)
                        for i = 0:na-1
                            u_range[i+1] = u_min + (i * (p1 - u_min)) / (na - 1)
                        end
                        u_range[na+1] = x.Pnl
                        nu = na + 1
                    else
                        u_range[1] = x.Pnl
                        nu = 1
                    end
                else
                    if p1 <= u_max && p1 >= u_min
                        na = min(Int(div(p1 - u_min, step) + 2), da)
                        for i = 0:na-1
                            u_range[i+1] = u_min + (i * (p1 - u_min)) / (na - 1)
                        end
                        nu = na
                    else
                        na = min(Int(div(u_max - u_min, step) + 2), da)
                        for i = 0:na-1
                            u_range[i+1] = u_min + (i * (u_max - u_min)) / (na - 1)
                        end
                        nu = na

                    end
                end

            end
        end
    end

    return nu
end


"""
    dynamic(x::State,Pbt::Float32,w::Float32,time_next::Vector{Float32},rep_days::Array{Float32, 2},rep_days_dev::Array{Float32, 2},mg)

TBW
"""
function dynamic(x, Pbt,Pnl_next, fc,el,hytank,bt,dt)

    # next net load submit to the uncertainty
    Pel_max, Pfc_max, Pel_min, Pfc_min = h2_bounds(x, fc,el,hytank,dt)
    #level of hydrogen update 
    Pnl = x.Pnl

    δ_next = 0
    LoH_next = x.LoH 
    PH2=0.0f0

    if (Pnl - Pbt > 1e-12)
        PH2 = min(Pfc_max, max(Pfc_min, Pnl - Pbt))
        LoH_next = x.LoH - PH2 * fc.consumption_slope * dt
        if (PH2>1e-12)
            δ_next = 1
        end

    elseif (Pnl - Pbt <= Pel_min + 1e-12)
        PH2 = max(Pel_max, Pnl - Pbt)
        LoH_next = x.LoH - PH2 / el.consumption_slope * dt
        
        if (PH2<1e-12)
            δ_next = -1
        end

    end

    #battery dynamic 
    E_bt_next = max(0.0f0, min(x.Ebatt - (Pbt + bt.loss_factor * abs(Pbt)) * dt, bt.energy_rated))
    x_next = State(Pnl_next, Float32(LoH_next), Float32(E_bt_next), Int32(δ_next))
    Pshed = Pnl - Pbt - PH2

    return x_next, Pshed,PH2

end

function m_pnl()

    t = [1.0, 1, 1]

    dt = mg.project.timestep
    # next time indexes
    Pnl_gen = zeros(Float32, 8760)
    w_law = Normal(0, 1)
    for i in eachindex(Pnl_gen)
        Pnl_gen[i] = rep_days[Int(t[3]), Int(t[1])+1] + rand(w_law) * rep_days_dev[Int(t[3]), Int(t[1])+1]
        if i != length(Pnl_gen)
            time_dyn_forward(t, dt, mi_l)
        end
    end

    return Pnl_gen

end

"""
    gen_mark_c(probs,cats,ini)

TBW
"""
function gen_mark_c(probs,cats,ini)

    Pnl=zeros(Float32,8760)
    Pnl[1]=ini
    for i=2:8760
        h=data.Hour[i]
        m=data.Month[i]
         h_prev=data.Hour[i-1]
        m_prev=data.Month[i-1]
       cat = Int(abs(div(Pnl[i-1] - cats[m_prev, h_prev+1, 1], cats[m_prev, h_prev+1, 2] - cats[m_prev, h_prev+1, 1])) + 1)

        if cat > 16
            cat = 16
        end
        prob_val=zeros(Float32,count(x->(x>0.0f0),probs[m_prev,h_prev+1,:,cat]))
        indexes=zeros(Int32,length(prob_val))
        z=0
        
        
            for j=1:length(probs[m_prev,h_prev+1,:,cat])
               
                if probs[m_prev,h_prev+1,j,cat] > 0.0f0
                    z=z+1
                   indexes[z]=j
                    prob_val[z]=probs[m_prev,h_prev+1,j,cat]
                end
            end
        Pnl[i]=cats[m,h+1,indexes[sample(Weights(prob_val))]+1]
    end
    return Pnl
end

"""
    CRF(i::Float32,T)

TBW
"""
function CRF(i::Float32, T)
    if i != 0.0f0
        a = (1 + i)^T
        return i * a / (a - 1)
    else
        return 1 / T
    end
end
CRF(0.05f0, 25.0f0)
CRFproj(T::Float32) = CRF(mg.project.discount_rate, T)


function dispatch_1(s::State, fc,el,hytank,bt,dt)

    Pbatt_cmax, Pbatt_dmax = bt_bounds(s, bt,dt)
    Pel_max, Pfc_max, Pel_min, Pfc_min = h2_bounds(s,fc,el,hytank,dt)
    Pnl, Pgen, Pbatt, Pspill, Pshed, Pdump, Pelyz, Pfc = Microgrids.dispatch_1(s.Pnl, Pbatt_cmax, Pbatt_dmax, 0.0f0, 0.0f0, Pel_min, Pel_max, Pfc_min, Pfc_max)

    return Pbatt, Pfc - Pelyz
end

function cost_to_go(pen::Float32,Pshed::Float32,Ph2::Float32,Pbt::Float32,δ::Int32,fc,el,bt,d_rt::Float32=deg_ratio_rt,d_st::Float32=deg_ratio_st,Eol::Float32=0.10f0)
    c=0
    
    c+= pen * max(0.0f0, Pshed)
    if Ph2> ( 10e-12)
            c += (d_rt/Eol)*(fc.power_rated*fc.investment_price)
        if δ !=1
            c += (d_st/Eol)*(fc.power_rated*fc.investment_price)
        end
    elseif (Ph2<  -10e-12)
         c += (d_rt/Eol)*(el.power_rated*el.investment_price)
        if δ !=-1
            c += (d_st/Eol)*(el.power_rated*el.investment_price)
        end
    end
    if Pbt<=-10e-12 || Pbt>=10e-12
        c+=(abs(Pbt)*bt.investment_price)/(2*bt.lifetime_cycles)
    end
    return c
end



"""
    DPrecursion(K::Int32,pen::Float32,time_init::Vector{Float32},max_pnl::Float32,min_pnl::Float32,mi_l::Vector{Int32},rep_days,rep_days_dev,mg)

TBW
"""
function DPrecursion(K::Int32, pen::Float32,  max_pnl::Float32, min_pnl::Float32, ref::Vector{Float32},σ::Float32,ϕ::Float32, mg,d_rt::Float32=deg_ratio_rt,d_st::Float32=deg_ratio_st)

   
    bt = mg.storage
    el = mg.electrolyzer[1]
    fc = mg.dispatchables.fuel_cell[1]
    hytank = mg.tanks.h2Tank
    dt = mg.project.timestep
    nPnl = 16
    nEbt = 8
    nloh = 8

    Pnl_grid = range(min_pnl, max_pnl, nPnl)
    Ebt_grid = range(bt.energy_rated * bt.SoC_min, bt.energy_rated * bt.SoC_max, nEbt)
    LoH_grid = range(hytank.capacity * hytank.min_filling_ratio, hytank.capacity * hytank.max_filling_ratio, nloh)
    δ_grid =    range(-1,1,3)

   
    IBS = BSpline(Linear())

    nw = 11
    w_law = Normal(0, σ)
    demand = range(-3σ, 3σ, nw)
    prob_demand = pdf.(w_law, demand) / (sum(pdf(w_law, demand)))

    Jopt = zeros(Float32, nPnl, nloh, nEbt,3, K + 1)
    # Policy over time
    Pbt_grid = zeros(nthreads(), 25)

    policy = zeros(Float16, nPnl, nloh, nEbt, 3, K)
    polh2 = zeros(Float16, nPnl, nloh, nEbt, 3, K)
    u_rb = zeros(Float16, nPnl, nloh, nEbt,3, 2)


    @inbounds for k = K:-1:1

        Jopt_next = Interpolations.scale(interpolate!(selectdim(Jopt, ndims(Jopt), k + 1), IBS), Pnl_grid, LoH_grid, Ebt_grid, δ_grid)

        @inbounds @threads for pnl = 1:nPnl

            @inbounds for loh = 1:nloh, ebt = 1:nEbt,δ =1:3

                xk = State(Pnl_grid[pnl], LoH_grid[loh], Ebt_grid[ebt], δ_grid[δ])

                nu = u_range!(selectdim(Pbt_grid, 1, threadid()), xk, fc,el,hytank,bt,dt, 25, 12, 12) #peut changer de taille à régler
                err=(xk.Pnl-ref[k])
                #CREATE UGRID
                Jk_xu_min = Inf
                u_opt = Pbt_grid[threadid(), 1]
                ph_opt = 0.0f0

                @inbounds for u_i = 1:nu  #gérer les situations où il n'y pas de u
                    uk = Pbt_grid[threadid(), u_i]
                    Jk_xu = 0.0f0

                    PH2 = 0.0f0
                     x_next, Pshed, PH2 = dynamic(xk, uk, 0.0f0,fc,el,hytank,bt,dt)

                    @inbounds @fastmath for iw in eachindex(demand) # for each uncertainty value                
                        Jk_xu += Jopt_next(ref[k+1] + (err * ϕ + demand[iw]) , x_next.LoH, x_next.Ebatt, x_next.δ) * prob_demand[iw]

                    end # each uncertainty

                     Jk_xu += cost_to_go(pen,Pshed,PH2,uk,xk.δ,fc,el,bt)  #+ max(x_next.δ,0)*(d_rt/1-Eol)*(fc.investment_price)
                    
                    
                    if (Jk_xu < Jk_xu_min)
                        u_opt = uk
                        Jk_xu_min = Jk_xu
                        ph_opt = PH2
                    end

                end #each control


                Jopt[pnl, loh, ebt, δ, k] = Jk_xu_min
                policy[pnl, loh, ebt, δ, k] = u_opt
                polh2[pnl, loh, ebt, δ, k] = ph_opt

            end # each state
        end
        #      
    end #each instant


    return Jopt, policy, polh2, u_rb

end


function DPrecursion_jopt(K::Int32, pen::Float32,ref::Vector{Float32},σ::Float32,ϕ::Float32, Pnl_cate,mg,d_rt::Float32=deg_ratio_rt,d_st::Float32=deg_ratio_st)

    el=mg.electrolyzer[1]
    fc= mg.dispatchables.fuel_cell[1]
    bt = mg.storage
    hytank = mg.tanks.h2Tank
    dt = mg.project.timestep
    nPnl = 40
    nEbt =40
    nloh = 40
    Pnl_grid = zeros(Float32, 12, 24, nPnl)
    for i = 1:12, j = 1:24
        Pnl_grid[i, j, :] = range(Pnl_cate[i, j, 1], Pnl_cate[i, j, end], nPnl)
    end
    Ebt_grid = range(bt.energy_rated * bt.SoC_min, bt.energy_rated * bt.SoC_max, nEbt)
    LoH_grid = range(hytank.capacity * hytank.min_filling_ratio, hytank.capacity * hytank.max_filling_ratio, nloh)
    δ_grid =    range(-1,1,3)

    time=[23.0,31,12]

    
    IBS = BSpline(Linear())

    nw = 11
    w_law = Normal(0, σ)
    demand = range(-3σ, 3σ, nw)
    prob_demand = pdf.(w_law, demand) / (sum(pdf(w_law, demand)))

    Jopt = zeros(Float32, nPnl, nloh, nEbt,3, K + 1)
    # Policy over time
    Pbt_grid = zeros(nthreads(), 50)


    @inbounds for k = K:-1:1
        time_next = time
        time_dyn_backward!(time, dt, mi_l)

        m = Int(time[3])
        h = Int(time[1])
        m_next = Int(time_next[3])
        h_next = Int(time_next[1])

        Jopt_next = Interpolations.scale(interpolate!(selectdim(Jopt, ndims(Jopt), k + 1), IBS),range(Pnl_cate[m_next, h_next+1, 1], Pnl_cate[m_next, h_next+1, end], nPnl), LoH_grid, Ebt_grid, δ_grid)

        @inbounds @threads for pnl = 1:nPnl

            @inbounds for loh = 1:nloh, ebt = 1:nEbt,δ =1:3

                xk = State(Pnl_grid[m, h+1, pnl], LoH_grid[loh], Ebt_grid[ebt], δ_grid[δ])

                nu = u_range!(selectdim(Pbt_grid, 1, threadid()), xk, fc,el,hytank,bt,dt, 25, 12, 12)  #peut changer de taille à régler
                err=(xk.Pnl-ref[k])
                #CREATE UGRID
                Jk_xu_min = Inf
                u_opt = Pbt_grid[threadid(), 1]


                @inbounds for u_i = 1:nu  #gérer les situations où il n'y pas de u
                    uk = Pbt_grid[threadid(), u_i]
                    Jk_xu = 0.0f0


                     x_next, Pshed, PH2 = dynamic(xk, uk, 0.0f0,fc,el,hytank,bt,dt)

                    @inbounds @fastmath for iw in eachindex(demand) # for each uncertainty value                
                        Jk_xu += Jopt_next(ref[k+1] + (err * ϕ + demand[iw]), x_next.LoH, x_next.Ebatt, x_next.δ) * prob_demand[iw]

                    end # each uncertainty

                     Jk_xu += cost_to_go(pen,Pshed,PH2,uk,xk.δ,fc,el,bt) 
                    
                    if (Jk_xu < Jk_xu_min)
                        u_opt = uk
                        Jk_xu_min = Jk_xu

                    end

                end #each control


                Jopt[pnl, loh, ebt, δ, k] = Jk_xu_min
            end # each state
        end
        #      

    end #each instant


    return Jopt

end
function DP_jopt(K::Int32, pen::Float32,Pnl::Vector{Float32},mg)

     bt = mg.storage
    el = mg.electrolyzer[1]
    fc = mg.dispatchables.fuel_cell[1]
    hytank = mg.tanks.h2Tank
    dt = mg.project.timestep
   
    nEbt =50
    nloh = 50
   
    Ebt_grid = range(bt.energy_rated * bt.SoC_min, bt.energy_rated * bt.SoC_max, nEbt)
    LoH_grid = range(hytank.capacity * hytank.min_filling_ratio, hytank.capacity * hytank.max_filling_ratio, nloh)
    δ_grid =    range(-1,1,3)

    

    
    IBS = BSpline(Linear())

   
    Jopt = zeros(Float32,nloh, nEbt,3, K + 1)
    # Policy over time
    Pbt_grid = zeros(nthreads(), 50)


    @inbounds for k = K:-1:1
       
       
        Jopt_next = Interpolations.scale(interpolate!(selectdim(Jopt, ndims(Jopt), k + 1), IBS), LoH_grid, Ebt_grid, δ_grid)

        @inbounds @threads for loh = 1:nloh

            @inbounds for ebt = 1:nEbt,δ =1:3

                xk = State(Pnl[k], LoH_grid[loh], Ebt_grid[ebt], δ_grid[δ])

                 nu = u_range!(selectdim(Pbt_grid, 1, threadid()), xk, fc,el,hytank,bt,dt, 25, 12, 12) #peut changer de taille à régler
                
                #CREATE UGRID
                Jk_xu_min = Inf
                u_opt = Pbt_grid[threadid(), 1]


                @inbounds for u_i = 1:nu  #gérer les situations où il n'y pas de u
                    uk = Pbt_grid[threadid(), u_i]
                    Jk_xu = 0.0f0


                     x_next, Pshed, PH2 = dynamic(xk, uk, 0.0f0,fc,el,hytank,bt,dt)

                               
                        Jk_xu = Jopt_next( x_next.LoH, x_next.Ebatt, x_next.δ) + cost_to_go(pen,Pshed,PH2,uk,xk.δ,fc,el,bt) 
                    
                    if (Jk_xu < Jk_xu_min)
                        u_opt = uk
                        Jk_xu_min = Jk_xu

                    end

                end #each control


                Jopt[loh, ebt, δ, k] = Jk_xu_min
            end # each state
        end
        #      

    end #each instant


    return Jopt

end

function DPrecursion_jopt_mc(K::Int32, pen::Float32, time_init::Vector{Float32}, mi_l::Vector{Int32}, Pnl_cate, Pnl_prob, mg;mode::Int32=Int32(1))

      bt = mg.storage
    el = mg.electrolyzer[1]
    fc = mg.dispatchables.fuel_cell[1]
    hytank = mg.tanks.h2Tank
    dt = mg.project.timestep
    nPnl = 40
    nEbt = 40
    nloh = 40
   
    

    n_cat = length(Pnl_cate[1, 1, :]) - 1
     Pnl_grid = zeros(Float32, 12, 24, nPnl)
    for i = 1:12, j = 1:24
        Pnl_grid[i, j, :] = range(Pnl_cate[i, j, 1], Pnl_cate[i, j, end], nPnl)
    end
    Ebt_grid = range(bt.energy_rated * bt.SoC_min, bt.energy_rated * bt.SoC_max, nEbt)
    LoH_grid = range(hytank.capacity * hytank.min_filling_ratio, hytank.capacity * hytank.max_filling_ratio, nloh)
    δ_grid =    range(-1,1,3)

    time = time_init

    IBS = BSpline(Linear())



    Jopt = zeros(Float32, nPnl, nloh, nEbt, 3, K + 1)
    # Policy over time
    Pbt_grid = zeros(Float32,nthreads(), 50)


    @inbounds for k = K:-1:1
        time_next = time
        time_dyn_backward!(time, dt, mi_l)

        m = Int(time[3])
        h = Int(time[1])
        m_next = Int(time_next[3])
        h_next = Int(time_next[1])
        Jopt_next = Interpolations.scale(interpolate!(selectdim(Jopt, ndims(Jopt), k + 1), IBS), range(Pnl_cate[m_next, h_next+1, 1], Pnl_cate[m_next, h_next+1, end], nPnl), LoH_grid, Ebt_grid, δ_grid)


        @inbounds for pnl = 1:nPnl

            @inbounds for loh = 1:nloh, ebt = 1:nEbt,δ=1:3

                xk = State(Pnl_grid[m, h+1, pnl], LoH_grid[loh], Ebt_grid[ebt],δ_grid[δ])
                cat = Int(div(xk.Pnl - Pnl_cate[m, h+1, 1], Pnl_cate[m, h+1, 2] - Pnl_cate[m, h+1, 1]) + 1)
                if cat > 16
                     cat = 16
                end
                 nu = u_range!(selectdim(Pbt_grid, 1, threadid()), xk, fc,el,hytank,bt,dt, 25, 12, 12)  #peut changer de taille à régler

                #CREATE UGRID
                Jk_xu_min = Inf32
                u_opt = Pbt_grid[threadid(), 1]


                @inbounds for u_i = 1:nu  #gérer les situations où il n'y pas de u
                    uk = Pbt_grid[threadid(), u_i]
                    Jk_xu = 0.0f0
                     x_next, Pshed, PH2 = dynamic(xk, uk, 0.0f0,fc,el,hytank,bt,dt)
                    
                   @inbounds @fastmath  for iw = 1:n_cat # for each uncertainty value  
                        if Pnl_prob[m, h+1, iw, cat] > 1e-10
                            if mode==1
                                Jk_xu += Jopt_next(Pnl_cate[m_next, h_next+1, iw] , x_next.LoH, x_next.Ebatt, x_next.δ) * Pnl_prob[m, h+1, iw, cat]
                            elseif mode==2
                                 Jk_xu += Jopt_next((Pnl_cate[m_next, h_next+1, iw] + Pnl_cate[m_next, h_next+1, iw+1])/2, x_next.LoH, x_next.Ebatt, x_next.δ) * Pnl_prob[m, h+1, iw, cat]
                            elseif mode==3
                                 Jk_xu += Jopt_next(Pnl_cate[m_next, h_next+1, iw+1] , x_next.LoH, x_next.Ebatt, x_next.δ) * Pnl_prob[m, h+1, iw, cat]
                            elseif mode==4
                                Jk_xu += (Jopt_next(Pnl_cate[m_next, h_next+1, iw+1] , x_next.LoH, x_next.Ebatt, x_next.δ) + Jopt_next(Pnl_cate[m_next, h_next+1, iw] , x_next.LoH, x_next.Ebatt, x_next.δ))  * Pnl_prob[m, h+1, iw, cat]/2
                            end
                                #  @assert Jk_xu >= 0.0f0 " problem, $m, $(h+1), $cat, $Jk_xu, $(Pnl_prob[m, h+1, iw, cat])"
                        end
                    end # each uncertainty

                      Jk_xu += cost_to_go(pen,Pshed,PH2,uk,xk.δ,fc,el,bt) 
                    if (Jk_xu <= Jk_xu_min)
                        u_opt = uk
                        Jk_xu_min = Jk_xu

                    end

                end #each control

              #  @assert Jk_xu_min >= Jopt[pnl, loh, ebt, δ, K+1] - 0.1 " problem, $m, $(h+1), $cat,$xk, $(dynamic(xk, u_opt, 0.0f0, mg)), $u_opt,$( Jk_xu_min), $(Jopt[pnl,loh,ebt,cy,rt1,rt2,K+1])"

                Jopt[pnl, loh, ebt, δ, k] = Jk_xu_min


            end # each state
        end
        #      


    end #each instant


    return Jopt

end

function DPrecursion_uopt(K::Int32, pen::Float32,  max_pnl::Float32, min_pnl::Float32, ref,σ,ϕ, mg)
   bt = mg.storage
    el = mg.electrolyzer[1]
    fc = mg.dispatchables.fuel_cell[1]
    hytank = mg.tanks.h2Tank
    dt = mg.project.timestep
    nPnl = 16
    nEbt = 8
    nloh = 8

    Pnl_grid = range(min_pnl, max_pnl, nPnl)
    Ebt_grid = range(bt.energy_rated * bt.SoC_min, bt.energy_rated * bt.SoC_max, nEbt)
    LoH_grid = range(hytank.capacity * hytank.min_filling_ratio, hytank.capacity * hytank.max_filling_ratio, nloh)
    δ_grid =    range(-1,1,3)



    
    IBS = BSpline(Linear())

    nw = 11
    w_law = Normal(0, σ)
    demand = range(-3σ, 3σ, nw)
    prob_demand = pdf.(w_law, demand) / (sum(pdf(w_law, demand)))

    Jopt = zeros(Float32, nPnl, nloh, nEbt,3, K + 1)
    # Policy over time
    Pbt_grid = zeros(nthreads(), 25)

    policy = zeros(Float16, nPnl, nloh, nEbt,3, K)

    u_rb = zeros(Float16, nPnl, nloh, nEbt,3, 2)

   
    @inbounds for k = K:-1:1

        Jopt_next = Interpolations.scale(interpolate!(selectdim(Jopt, ndims(Jopt), k + 1), IBS), Pnl_grid, LoH_grid, Ebt_grid, δ_grid)

        @inbounds @threads for pnl = 1:nPnl

            @inbounds for loh = 1:nloh, ebt = 1:nEbt, δ=1:3

                xk = State(Pnl_grid[pnl], LoH_grid[loh], Ebt_grid[ebt], δ_grid[δ])
                 err=(xk.Pnl-ref[k])
                 nu = u_range!(selectdim(Pbt_grid, 1, threadid()), xk, fc,el,hytank,bt,dt, 25, 12, 12) #peut changer de taille à régler

                #CREATE UGRID
                Jk_xu_min = Inf
                u_opt = Pbt_grid[threadid(), 1]


                @inbounds for u_i = 1:nu  #gérer les situations où il n'y pas de u
                    uk = Pbt_grid[threadid(), u_i]
                    Jk_xu = 0.0f0


                      x_next, Pshed, PH2 = dynamic(xk, uk, 0.0f0,fc,el,hytank,bt,dt)

                    @inbounds @fastmath for iw in eachindex(demand) # for each uncertainty value                
                        Jk_xu += Jopt_next(ref[k+1] + (err * ϕ + demand[iw]) , x_next.LoH, x_next.Ebatt, x_next.δ) * prob_demand[iw]

                    end # each uncertainty

                    Jk_xu += cost_to_go(pen,Pshed,PH2,uk,xk.δ,fc,el,bt) 
                    if (Jk_xu < Jk_xu_min)
                        u_opt = uk
                        Jk_xu_min = Jk_xu

                    end

                end #each control


                Jopt[pnl, loh, ebt, δ, k] = Jk_xu_min
                policy[pnl, loh, ebt, δ, k] = u_opt

            end # each state
        end
        #      

       

    end #each instant


    return policy, u_rb

end


function DPrecursion_rb_cost(K::Int32, pen::Float32, max_pnl::Float32, min_pnl::Float32,ref,σ,ϕ, mg)
    
      bt = mg.storage
    el = mg.electrolyzer[1]
    fc = mg.dispatchables.fuel_cell[1]
    hytank = mg.tanks.h2Tank
    dt = mg.project.timestep
    nPnl = 16
    nEbt = 8
    nloh = 8

    Pnl_grid = range(min_pnl, max_pnl, nPnl)
    Ebt_grid = range(bt.energy_rated * bt.SoC_min, bt.energy_rated * bt.SoC_max, nEbt)
    LoH_grid = range(hytank.capacity * hytank.min_filling_ratio, hytank.capacity * hytank.max_filling_ratio, nloh)
    δ_grid =    range(-1,1,3)



    IBS = BSpline(Linear())

    nw = 11
     w_law = Normal(0, σ)
    demand = range(-3σ, 3σ, nw)
    prob_demand = pdf.(w_law, demand) / (sum(pdf(w_law, demand)))

    Jopt = zeros(Float32, nPnl, nloh, nEbt, 3, K + 1)
    # Policy over time




    u_rb = zeros(Float16, nPnl, nloh, nEbt, 3, 2)



    @inbounds for k = K:-1:1

         Jopt_next = Interpolations.scale(interpolate!(selectdim(Jopt, ndims(Jopt), k + 1), IBS), Pnl_grid, LoH_grid, Ebt_grid, δ_grid)

        @inbounds @threads for pnl = 1:nPnl

            @inbounds for loh = 1:nloh, ebt = 1:nEbt, δ=1:3

                xk = State(Pnl_grid[pnl], LoH_grid[loh], Ebt_grid[ebt], δ_grid[δ])
                 err=(xk.Pnl-ref[k])

                Jk_xu = 0.0f0
                uk = u_rb[pnl, loh, ebt, δ, 1]

                 x_next, Pshed, PH2 = dynamic(xk, uk, 0.0f0,fc,el,hytank,bt,dt)

                @inbounds @fastmath for iw in eachindex(demand) # for each uncertainty value                
                Jk_xu += Jopt_next(ref[k+1] + (err * ϕ + demand[iw]) , x_next.LoH, x_next.Ebatt, x_next.δ) * prob_demand[iw]
                end # each uncertainty

                Jk_xu += cost_to_go(pen,Pshed,PH2,uk,xk.δ,fc,el,bt)  


                Jopt[pnl, loh, ebt, δ, k] = Jk_xu


            end # each state
        end
        #      

    end #each instant


    return Jopt[:, :, :, :, 1]

end
function sdp_simu_jopt(Jopt, n, pen,ϕ,σ,Pnl_cate, mg1,ref::Vector{Float32},Pnl::Vector{Float32},data=data)
   
  
     bt = mg1.storage
    el = mg1.electrolyzer[1]
    fc = mg1.dispatchables.fuel_cell[1]
    hytank = mg1.tanks.h2Tank
    dt = mg1.project.timestep
    Pren = zeros(n)
    Pren = production(mg1.nondispatchables[1]) .+ production(mg1.nondispatchables[2])
    Pbt = zeros(n)
    Pel = zeros(n)
    Pfc = zeros(n)
    Ebt = zeros(n + 1)
    Loh = zeros(n + 1)
    Lof = zeros(n + 1)
    Pgen = zeros(n)
    Phb = zeros(n)
    Pshed = zeros(n)
    Pdump = zeros(n)
    Pspill = zeros(n)
   
    IBS = BSpline(Linear())
    δ=0
    Ebt[1] = 0.0f0
    Loh[1] = hytank.ini_filling_ratio * hytank.capacity
    nEbt = size(Jopt)[3]
    nPnl = size(Jopt)[1]
    nloh = size(Jopt)[2]
    
    Ebt_grid = range(bt.energy_rated * bt.SoC_min, bt.energy_rated * bt.SoC_max, nEbt)
    Pnl_grid = zeros(Float32, 12, 24, nPnl)

    for i = 1:12, j = 1:24
        Pnl_grid[i, j, :] = range(Pnl_cate[i, j, 1], Pnl_cate[i, j, end], nPnl)
    end

    LoH_grid = range(hytank.capacity * hytank.min_filling_ratio, hytank.capacity * hytank.max_filling_ratio, nloh)
    δ_grid = range(-1,1,3)
    
    Pbt_grid = zeros(100)
    nw = 11
    w_law = Normal(0, σ)
    demand = range(-3σ, 3σ, nw)
    prob_demand = pdf.(w_law, demand) / (sum(pdf(w_law, demand)))
    for k = 1:n
         if k < n
            time_next = [Float32(data.Hour[k+1]), Float32(data.Day[k+1]), Float32(data.Month[k+1])]

        else
            time_next = [0.0f0, 1, 1]
        end

        time = [Float32(data.Hour[k]), Float32(data.Day[k]), Float32(data.Month[k])]
        m = Int(time[3])
        h = Int(time[1])
        m_next = Int(time_next[3])
        h_next = Int(time_next[1])
        J_interp = Interpolations.extrapolate(Interpolations.scale(Interpolations.interpolate!(selectdim(Jopt, ndims(Jopt), k + 1), IBS), range(Pnl_cate[m_next, h_next+1, 1], Pnl_cate[m_next, h_next+1, end], nPnl), LoH_grid, Ebt_grid, δ_grid), Line())

        xk = State(Pnl[k], Loh[k], Ebt[k], δ)
        err=(xk.Pnl-ref[k])
        nu = u_range!(Pbt_grid, xk,fc,el,hytank,bt,dt, 5, 49, 49)
       
        
        Jk_xu_min = Inf
        u_opt = 0.0f0
        ph_opt = 0.0f0

        @inbounds for u_i = 1:nu  #gérer les situations où il n'y pas de u
            uk = Pbt_grid[u_i]
            Jk_xu = 0.0f0

           x_next, Ps, Ph2 = dynamic(xk, uk, 0.0f0, fc,el,hytank,bt,dt)

                @inbounds @fastmath for iw in eachindex(demand) # for each uncertainty value                
                Jk_xu += J_interp(ref[k+1] + (err * ϕ + demand[iw]) , x_next.LoH, x_next.Ebatt, x_next.δ) * prob_demand[iw]
                end # each uncertainty



             Jk_xu += cost_to_go(pen,Ps,Ph2,uk,xk.δ,fc,el,bt)
            if (Jk_xu < Jk_xu_min)
                u_opt = uk
                Jk_xu_min = Jk_xu
                ph_opt = Ph2
            end

        end #each control
        x_next, Ps, Ph2 = dynamic(xk, uk, 0.0f0, fc,el,hytank,bt,dt)
        Pbt[k] = u_opt
        Pel[k] = ifelse(Ph2 > 0.0f0, 0.0f0, -Ph2)
        Pfc[k] = ifelse(Ph2 > 0.0f0, Ph2, 0.0f0)
        Ebt[k+1] = x_next.Ebatt
        Loh[k+1] = x_next.LoH
        Pspill[k] = ifelse(Ps > 0.0f0, 0.0f0, -Ps)
        Pshed[k] = ifelse(Ps > 0.0f0, Ps, 0.0f0)
         δ=x_next.δ

    end
    traj = OperationTraj(Pnl, Pshed, Pren, Pgen, Pfc, Pel, Phb, Ebt, Pbt, Loh, Lof, Pspill, Pdump)
    stats = aggregation(mg1, traj)
    # Eval the microgrid costs
    costs = economics(mg1, stats)

    return traj, stats, costs

end

function Dp_simu(Jopt, n, pen, mg1,Pnl::Vector{Float32},data=data)
   
   
     bt = mg1.storage
    el = mg1.electrolyzer[1]
    fc = mg1.dispatchables.fuel_cell[1]
    hytank = mg1.tanks.h2Tank
    dt = mg1.project.timestep
    Pren = zeros(n)
    Pren = production(mg1.nondispatchables[1]) .+ production(mg1.nondispatchables[2])
    Pbt = zeros(n)
    Pel = zeros(n)
    Pfc = zeros(n)
    Ebt = zeros(n + 1)
    Loh = zeros(n + 1)
    Lof = zeros(n + 1)
    Pgen = zeros(n)
    Phb = zeros(n)
    Pshed = zeros(n)
    Pdump = zeros(n)
    Pspill = zeros(n)
   
    IBS = BSpline(Linear())
    δ=0
    Ebt[1] = 0.0f0
    Loh[1] = hytank.ini_filling_ratio * hytank.capacity
    nEbt = size(Jopt)[2]
    
    nloh = size(Jopt)[1]
    
    Ebt_grid = range(bt.energy_rated * bt.SoC_min, bt.energy_rated * bt.SoC_max, nEbt)
  

   
    LoH_grid = range(hytank.capacity * hytank.min_filling_ratio, hytank.capacity * hytank.max_filling_ratio, nloh)
    δ_grid = range(-1,1,3)
    
    Pbt_grid = zeros(100)
   
    for k = 1:n
         if k < n
            time_next = [Float32(data.Hour[k+1]), Float32(data.Day[k+1]), Float32(data.Month[k+1])]

        else
            time_next = [0.0f0, 1, 1]
        end

        time = [Float32(data.Hour[k]), Float32(data.Day[k]), Float32(data.Month[k])]
        m = Int(time[3])
        h = Int(time[1])
        m_next = Int(time_next[3])
        h_next = Int(time_next[1])
        J_interp = Interpolations.extrapolate(Interpolations.scale(Interpolations.interpolate!(selectdim(Jopt, ndims(Jopt), k + 1), IBS), LoH_grid, Ebt_grid, δ_grid), Line())

        xk = State(Pnl[k], Loh[k], Ebt[k], δ)
      
         nu = u_range!(Pbt_grid, xk,fc,el,hytank,bt,dt, 5, 49, 49)
       
        
        Jk_xu_min = Inf
        u_opt = 0.0f0
        ph_opt = 0.0f0

        @inbounds for u_i = 1:nu  #gérer les situations où il n'y pas de u
            uk = Pbt_grid[u_i]
            Jk_xu = 0.0f0

           x_next, Ps, Ph2 = dynamic(xk, uk, 0.0f0, fc,el,hytank,bt,dt)

                          
                Jk_xu = J_interp( x_next.LoH, x_next.Ebatt, x_next.δ) + cost_to_go(pen,Ps,Ph2,uk,xk.δ,fc,el,bt)
                
            if (Jk_xu < Jk_xu_min)
                u_opt = uk
                Jk_xu_min = Jk_xu
                ph_opt = Ph2
            end

        end #each control
       x_next, Ps, Ph2 = dynamic(xk, uk, 0.0f0, fc,el,hytank,bt,dt)
        Pbt[k] = u_opt
        Pel[k] = ifelse(Ph2 > 0.0f0, 0.0f0, -Ph2)
        Pfc[k] = ifelse(Ph2 > 0.0f0, Ph2, 0.0f0)
        Ebt[k+1] = x_next.Ebatt
        Loh[k+1] = x_next.LoH
        Pspill[k] = ifelse(Ps > 0.0f0, 0.0f0, -Ps)
        Pshed[k] = ifelse(Ps > 0.0f0, Ps, 0.0f0)
         δ=x_next.δ

    end
    traj = OperationTraj(Pnl, Pshed, Pren, Pgen, Pfc, Pel, Phb, Ebt, Pbt, Loh, Lof, Pspill, Pdump)
    stats = aggregation(mg1, traj)
    # Eval the microgrid costs
    costs = economics(mg1, stats)

    return traj, stats, costs

end


function sdp_simu_jopt_mc(Jopt, n, pen, mg1, Pnl_cate, Pnl_prob,Pnl;mode::Int32=Int32(1),data=data)
   
   
     bt = mg1.storage
    el = mg1.electrolyzer[1]
    fc = mg1.dispatchables.fuel_cell[1]
    hytank = mg1.tanks.h2Tank
    dt = mg1.project.timestep
    Pren = zeros(n)
   Pren = production(mg1.nondispatchables[1]) .+ production(mg1.nondispatchables[2])

    Pbt = zeros(n)
    Pel = zeros(n)
    Pfc = zeros(n)
    Ebt = zeros(n + 1)
   
    Loh = zeros(n + 1)
    Lof = zeros(n + 1)
    Pgen = zeros(n)
    Phb = zeros(n)
    Pshed = zeros(n)
    Pdump = zeros(n)
    Pspill = zeros(n)
    
    IBS = BSpline(Linear())
    δ=0
    Ebt[1] = 0.0f0
    Loh[1] = hytank.ini_filling_ratio * hytank.capacity
    nEbt = size(Jopt)[3]
    nPnl = size(Jopt)[1]
    nloh = size(Jopt)[2]
   
    Ebt_grid = range(bt.energy_rated * bt.SoC_min, bt.energy_rated * bt.SoC_max, nEbt)
    Pnl_grid = zeros(Float32, 12, 24, nPnl)
    n_cat = length(Pnl_cate[1, 1, :]) - 1


    for i = 1:12, j = 1:24
        Pnl_grid[i, j, :] = range(Pnl_cate[i, j, 1], Pnl_cate[i, j, end], nPnl)
    end

    LoH_grid = range(hytank.capacity * hytank.min_filling_ratio, hytank.capacity * hytank.max_filling_ratio, nloh)
    δ_grid = range(-1,1,3)
    Pbt_grid = zeros(100)

    for k = 1:n
        if k < n
            time_next = [Float32(data.Hour[k+1]), Float32(data.Day[k+1]), Float32(data.Month[k+1])]

        else
            time_next = [0.0f0, 1, 1]
        end

        time = [Float32(data.Hour[k]), Float32(data.Day[k]), Float32(data.Month[k])]
        m = Int(time[3])
        h = Int(time[1])
        m_next = Int(time_next[3])
        h_next = Int(time_next[1])


        J_interp = Interpolations.extrapolate(Interpolations.scale(Interpolations.interpolate!(selectdim(Jopt, ndims(Jopt), k + 1), IBS), range(Pnl_cate[m_next, h_next+1, 1], Pnl_cate[m_next, h_next+1, end], nPnl),
                LoH_grid, Ebt_grid, δ_grid), Line())

        xk = State(Pnl[k], Float32(Loh[k]), Float32(Ebt[k]),Int32(δ))
        cat = Int(div(xk.Pnl - Pnl_cate[m, h+1, 1], Pnl_cate[m, h+1, 2] - Pnl_cate[m, h+1, 1]) + 1)

        if cat > 16
            cat = 16
        end
        
        if cat < 1
            cat = 1
        end
        
      #  @assert cat >= 1 && cat <= 16 "error $cat, $(xk.Pnl),  $(h+1), $m, $(Pnl_cate[m,h+1,1]), $(xk.Pnl-Pnl_cate[m,h+1,1]), $(Pnl_cate[m,h+1,2]-Pnl_cate[m,h+1,1]))"
         nu = u_range!(Pbt_grid, xk,fc,el,hytank,bt,dt, 5, 49, 49)
        #  nu_bt= Int(ceil((umax.P_bt-umin.P_bt)/dbt)) + 1
        #  nu_h2= Int(ceil((umax.P_H2-umin.P_H2)/dh2)) + 1




        #CREATE UGRID
        Jk_xu_min = Inf
        u_opt = 0.0f0
        ph_opt = 0.0f0


        @inbounds for u_i = 1:nu  #gérer les situations où il n'y pas de u
            uk = Pbt_grid[u_i]
            Jk_xu = 0.0f0
            x_next, Ps, Ph2 = dynamic(xk, uk, 0.0f0, fc,el,hytank,bt,dt)

            @fastmath for iw = 1:n_cat # for each uncertainty value  
                if Pnl_prob[m, h+1, iw, cat] != 0.0f0
                    if mode==1
                        Jk_xu += J_interp(Pnl_cate[m_next, h_next+1, iw], x_next.LoH, x_next.Ebatt, x_next.δ) * Pnl_prob[m, h+1, iw, cat]
                    elseif mode==2
                        Jk_xu += J_interp((Pnl_cate[m_next, h_next+1, iw]+Pnl_cate[m_next, h_next+1, iw+1])/2, x_next.LoH, x_next.Ebatt, x_next.δ) * Pnl_prob[m, h+1, iw, cat]
                    elseif mode==3 
                         Jk_xu += J_interp(Pnl_cate[m_next, h_next+1, iw+1], x_next.LoH, x_next.Ebatt, x_next.δ) * Pnl_prob[m, h+1, iw, cat]
                    elseif mode==4
                        Jk_xu += (J_interp(Pnl_cate[m_next, h_next+1, iw], x_next.LoH, x_next.Ebatt, x_next.δ)+J_interp(Pnl_cate[m_next, h_next+1, iw+1], x_next.LoH, x_next.Ebatt, x_next.δ)) * Pnl_prob[m, h+1, iw, cat]/2
                    end
                end
            end # each uncertainty

            #   @assert Jk_xu=0.0f0!=0 " cout nul"
          Jk_xu += cost_to_go(pen,Float32(Ps),Float32(Ph2),Float32(uk),xk.δ,fc,el,bt)
            if (Jk_xu <= Jk_xu_min)
                u_opt = uk
                Jk_xu_min = Jk_xu
                ph_opt = Ph2
            end

        end #each control
        x_next, Ps, Ph2 = dynamic(xk, u_opt, 0.0f0, fc,el,hytank,bt,dt)
        Pbt[k] = u_opt
        Pel[k] = ifelse(Ph2 > 0.0f0, 0.0f0, -Ph2)
        Pfc[k] = ifelse(Ph2 > 0.0f0, Ph2, 0.0f0)
        Ebt[k+1] = x_next.Ebatt
        Loh[k+1] = x_next.LoH
        Pspill[k] = ifelse(Ps > 0.0f0, 0.0f0, -Ps)
        Pshed[k] = ifelse(Ps > 0.0f0, Ps, 0.0f0)
        δ=x_next.δ


    end
    traj = OperationTraj(Float32.(Pnl), Float32.(Pshed), Float32.(Pren), Float32.(Pgen), Float32.(Pfc), Float32.(Pel), Float32.(Phb), Float32.(Ebt), Float32.(Pbt), Float32.(Loh), Float32.(Lof), Float32.(Pspill), Float32.(Pdump))
    stats = aggregation(mg1, traj)
    # Eval the microgrid costs
    costs = economics(mg1, stats)

    return traj, stats, costs
end


function moving_average(vs,q;even=true)

   
    n=length(vs)
    mav=zeros(n)
     vs_im=zeros(n+2q)
    vs_im[1:q].=vs[end-q+1:end]
    vs_im[n+q+1:n+2q].=vs[1:q]
    vs_im[q+1:n+q].=vs
    coef=ones(2q+1)
    if even
    coef[1]=0.5
    coef[end]=0.5
    end
        for i=1:n
            mav[i]=sum( vs_im[i:i+2q] .* coef) /sum(coef)
        end

    return mav
end

function sp_tr_remove(signal,per,d)
    n=length(signal)
    tre=zeros(length(signal))
    sig=zeros(n+2d*per)
    sig[1:d*per].=signal[end-d*per+1:end]
    sig[d*per+1:n+d*per] .= signal[:]
    sig[n+d*per+1:end]=signal[1:d*per]

    stv=zeros(n)
    vec=zeros(2d+1)
    for i=1:n
        tot=signal[i]
        for j=1:d
            tot+=sig[i-j*per+d*per] + sig[i+j*per+d*per] 
            vec[2(j-1)+1]=sig[i-j*per+d*per]
            vec[2j]=sig[i+j*per+d*per] 
        end
        vec[end]=signal[i]
        tre[i]=tot/(2d+1)
        stv[i]=std(vec)
    end
    return tre,stv

end

function hma_spe(signal,n)
    m=Int32((n-1)/2)
    new_sig=zeros(length(signal)+2m)
    new_sig[1:m]=signal[end-m+1:end]
    new_sig[m+1:end-m]=signal[:]
    new_sig[end-m+1:end]=signal[1:m]

    return hma(new_sig,n)[m+1:end-m]

end
function plot_periodigram(signal)
    n=length(signal)
    w=range(1,n/2,n)./n
    i=zeros(Float32,n)
    t=range(1,n)
    for k=1:n
        i[k]=(n^-1)*abs(sum(signal .* exp.(-2 .*pi .* im .*(t.-1) .*w[k])))^2
            
    end
   
    plot(w,i)
   # xscale("log")
    return i
end
function cal_w(a,j)
    return (1-a^(2j))*(1+a^2)*((1-a^2)^-1)-2j*a^(2j)
end
function season_estim(signal,d,q)
    
    n=length(signal)
    w=zeros(d)
    w_final=zeros(d)
    nj=Int32(n/d)
    for i=1:d
        if i<=q
            j=range(1,nj-1)
        else
             j=range(0,nj-2)
        end

        sum_int=0
        for z in j
            sum_int+=signal[i+z*d]
        end
        w[i]=sum_int/length(j)

    end
    for i in eachindex(w)
        w_final[i]=w[i]-mean(w)
    end
    return w_final
end