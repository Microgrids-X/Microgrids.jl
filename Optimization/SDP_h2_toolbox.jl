# Load Microgrid prior using this file
using Interpolations, Base.Threads, LinearAlgebra, Base.Iterators, Distributions,FreqTables,StatsBase
using Microgrids

# Structure describing the state of the Microgrid
struct State
    Pnl::Float64
    LoH::Float64
    Ebatt::Float64
    cycles::Float64
    RT_fc::Float64
    RT_el::Float64
end
"""
    State(a::NTuple{6, Float64})

TBW
"""
function State(a::NTuple{6,Float64})
    State(a[1], a[2], a[3], a[4], a[5], a[6])
end

struct Time
    hi::Float64
    di::Int64
    mi::Int64
end

"""
    time_dyn_forward(time::Array{Float64, 1},dt::Float64,mi_l::Vector{Int64})


"""
function time_dyn_forward(time::Array{Float64,1}, dt::Float64, mi_l::Vector{Int64})

    if time[1] < 23.0
        next_hi = time[1] + dt
        next_di = time[2]
        next_mi = time[3]
    else
        next_hi = 0.0
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
    time_dyn_backward!(time::Array{Float64, 1},dt::Float64,mi_l::Vector{Int64})

TBW
"""
function time_dyn_backward!(time::Array{Float64,1}, dt::Float64, mi_l::Vector{Int64})

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
function bt_bounds(x::State, mg::Microgrid)
    bt = mg.storage
    dt = mg.project.timestep
    if x.cycles >= mg.storage.lifetime_cycles / 5.0
        Pdis_max = 0.
        Pcha_max = 0.
    else
        Pdis_max = min((x.Ebatt - bt.energy_rated * bt.SoC_min) / ((1 + bt.loss_factor) * dt), bt.energy_rated * bt.discharge_rate)
        Pcha_max = max((x.Ebatt - bt.energy_rated * bt.SoC_max) / ((1 - bt.loss_factor) * dt), -bt.energy_rated * bt.charge_rate)
    end
    return Pcha_max, Pdis_max
end



"""
    h2_bounds( x::State,mg::Microgrid)

TBW
"""
function h2_bounds(x::State, mg::Microgrid)

    el = mg.electrolyzer[1]
    fc = mg.dispatchables.fuel_cell[1]
    hytank = mg.tanks.h2Tank
    dt = mg.project.timestep
    Pel_min, Pfc_min = -el.power_rated * el.minimum_load_ratio, fc.power_rated * fc.minimum_load_ratio

    if x.RT_el == el.lifetime_hours / 5.0

        Pel_max = 0.
    else
        Pel_max = ifelse((x.LoH - hytank.capacity * hytank.max_filling_ratio) * el.consumption_slope / dt > Pel_min, 0.0, max((x.LoH - hytank.capacity * hytank.max_filling_ratio) * el.consumption_slope / dt, -el.power_rated))
    end

    if x.RT_fc == fc.lifetime_hours / 5.0
        Pfc_max = 0.
    else
        Pfc_max = ifelse((x.LoH - hytank.capacity * hytank.min_filling_ratio) / (fc.consumption_slope * dt) < Pfc_min, 0.0, min((x.LoH - hytank.capacity * hytank.min_filling_ratio) / (fc.consumption_slope * dt), fc.power_rated))
    end

    return Pel_max, Pfc_max, Pel_min, Pfc_min
end

"""
    u_bounds( x::State,mg)

TBW
"""
function u_bounds(x::State, mg)
    bt = mg.storage
    dt = mg.project.timestep
    Pel_max, Pfc_max, Pel_min, Pfc_min = h2_bounds(x, mg)
    Pcha_max, Pdis_max = bt_bounds(x, mg)

    if x.Pnl >= 0.
        u_max = min(Pdis_max, x.Pnl)
        u_min = ifelse(Pfc_max >= x.Pnl, max(Pcha_max, x.Pnl - Pfc_max), min(Pdis_max, x.Pnl - Pfc_max))
    else
        u_max = max(min(x.Pnl - Pel_max, 0.0), Pcha_max)
        u_min = max(Pcha_max, x.Pnl - Pfc_max)
    end

    return u_min, u_max, Pel_min, Pfc_min
end


"""
    u_range!(u_range,x::State,mg,step=100,db=10,da=10)

TBW
"""
function u_range!(u_range, x::State, mg, step=100, db=10, da=10)

    u_min, u_max, Pel_min, Pfc_min = u_bounds(x, mg)

    na = da
    nb = db
    p1 = x.Pnl - Pfc_min
    p2 = x.Pnl - Pel_min
    nu = 0
    #  u_range=zeros(25)

    if x.Pnl >= 0.
        if u_min != u_max

            if p1 <= u_max && p1 >= u_min
                if p1 >= 0.

                    if 0. <= u_max && 0. >= u_min
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
                    if 0. <= u_max && 0. >= u_min
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
    dynamic(x::State,Pbt::Float64,w::Float64,time_next::Vector{Float64},rep_days::Array{Float64, 2},rep_days_dev::Array{Float64, 2},mg)

TBW
"""
function dynamic(x::State, Pbt::Float64,Pnl_next::Float64, mg)

    bt = mg.storage

    el = mg.electrolyzer[1]
    fc = mg.dispatchables.fuel_cell[1]

    dt = mg.project.timestep
    # next net load submit to the uncertainty
    Pel_max, Pfc_max, Pel_min, Pfc_min = h2_bounds(x, mg)
    #level of hydrogen update 
    Pnl = x.Pnl

    next_RT_el = x.RT_el
    next_RT_fc = x.RT_fc

    if (Pnl - Pbt > 1e-12)
        PH2 = min(Pfc_max, max(Pfc_min, Pnl - Pbt))
        LoH_next = x.LoH - PH2 * fc.consumption_slope * dt
        if (PH2>1e-12)
            next_RT_fc = x.RT_fc + dt
        end

    elseif (Pnl - Pbt <= Pel_min + 1e-12)
        PH2 = max(Pel_max, Pnl - Pbt)
        LoH_next = x.LoH - PH2 / el.consumption_slope * dt
        
        next_RT_el = x.RT_el + dt

    else
        PH2 = 0.
        LoH_next = x.LoH
    end

    #battery dynamic 
    E_bt_next = max(0.0, min(x.Ebatt - (Pbt + bt.loss_factor * abs(Pbt)) * dt, bt.energy_rated))
    cycles_next = x.cycles + (abs(Pbt) * dt) / (2 * bt.lifetime_cycles)

    x_next = State(Pnl_next, LoH_next, E_bt_next, cycles_next, next_RT_fc, next_RT_el)
    Pshed = Pnl - Pbt - PH2

    return x_next, Pshed,PH2

end

function m_pnl()

    t = [1.0, 1, 1]

    dt = mg.project.timestep
    # next time indexes
    Pnl_gen = zeros(Float64, 8760)
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

    Pnl=zeros(Float64,8760)
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
        prob_val=zeros(Float64,count(x->(x>0.),probs[m_prev,h_prev+1,:,cat]))
        indexes=zeros(Int64,length(prob_val))
        z=0
        
        
            for j=1:length(probs[m_prev,h_prev+1,:,cat])
               
                if probs[m_prev,h_prev+1,j,cat] > 0.
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
    simple_dyn(x::State,Pbt::Float64,mg)

TBW
"""
function simple_dyn(x::State, Pbt::Float64, mg)

    bt = mg.storage
    el = mg.electrolyzer[1]
    fc = mg.dispatchables.fuel_cell[1]

    dt = mg.project.timestep

    Pel_max, Pfc_max, Pel_min, Pfc_min = h2_bounds(x, mg)

    Pnl = x.Pnl

    next_RT_el = x.RT_el
    next_RT_fc = x.RT_fc

    if (Pnl - Pbt > 0.)
        PH2 = min(Pfc_max, max(Pfc_min, Pnl - Pbt))
        LoH_next = x.LoH - PH2 * fc.consumption_slope * dt

        next_RT_fc = x.RT_fc + dt

    elseif (Pnl - Pbt <= Pel_min + 1e-12)
        PH2 = max(Pel_max, Pnl - Pbt)
        LoH_next = x.LoH - PH2 / el.consumption_slope * dt

        next_RT_el = x.RT_el + dt

    else
        PH2 = 0.
        LoH_next = x.LoH
    end

    #battery dynamic 
    E_bt_next = max(0.0, min(x.Ebatt - (Pbt + bt.loss_factor * abs(Pbt)) * dt, bt.energy_rated))
    cycles_next = x.cycles + (abs(Pbt) * dt) / (2 * bt.lifetime_cycles)

    x_next = State(0.0, LoH_next, E_bt_next, cycles_next, next_RT_fc, next_RT_el)
    Ps = Pnl - Pbt - PH2

    return x_next, Ps, PH2

end

"""
    CRF(i::Float64,T)

TBW
"""
function CRF(i::Float64, T)
    if i != 0.0
        a = (1 + i)^T
        return i * a / (a - 1)
    else
        return 1 / T
    end
end
CRF(0.05, 25.0)
CRFproj(T::Float64) = CRF(mg.project.discount_rate, T)

"""
    term(x::State,T,mg::Microgrid)

TBW
"""
function term(x::State, T, mg::Microgrid)

    bt = mg.storage
    el = mg.electrolyzer[1]
    fc = mg.dispatchables.fuel_cell[1]
    τ = mg.project.discount_rate
    hprice=mg.tanks.h2Tank.combustible_price
    co = T / 8760
    Lt_bt = min(bt.lifetime_calendar, co * (bt.lifetime_cycles / x.cycles))
    Lt_fc = min(fc.lifetime_calendar, co * (fc.lifetime_hours / x.RT_fc))
    Lt_el = min(el.lifetime_calendar, co * (el.lifetime_hours / x.RT_el))

    c = bt.investment_price * bt.energy_rated * CRF(τ, Lt_bt) +
        fc.investment_price * fc.power_rated * CRF(τ, Lt_fc) +
        el.investment_price * el.power_rated * CRF(τ, Lt_el) #+ hprice * x.LoH

    return c


end


function dispatch_1(s::State, mg::Microgrid)

    Pbatt_cmax, Pbatt_dmax = bt_bounds(s, mg)
    Pel_max, Pfc_max, Pel_min, Pfc_min = h2_bounds(s, mg)
    Pnl, Pgen, Pbatt, Pspill, Pshed, Pdump, Pelyz, Pfc = Microgrids.dispatch_1(s.Pnl, Pbatt_cmax, Pbatt_dmax, 0., 0., Pel_min, Pel_max, Pfc_min, Pfc_max)

    return Pbatt, Pfc - Pelyz
end


"""
    DPrecursion(K::Int64,pen::Float64,time_init::Vector{Float64},max_pnl::Float64,min_pnl::Float64,mi_l::Vector{Int64},rep_days,rep_days_dev,mg)

TBW
"""
function DPrecursion(K::Int64, pen::Float64,  max_pnl::Float64, min_pnl::Float64, ref,σ,ϕ, mg)

    bt = mg.storage
    # el=mg.electrolyzer[1]
    #fc= mg.dispatchables.fuel_cell[1]
    hytank = mg.tanks.h2Tank
    dt = mg.project.timestep
    nPnl = 16
    nEbt = 8
    nloh = 8
    nRT = 10
    ncycles = 10

    Pnl_grid = range(min_pnl, max_pnl, nPnl)
    Ebt_grid = range(bt.energy_rated * bt.SoC_min, bt.energy_rated * bt.SoC_max, nEbt)
    LoH_grid = range(hytank.capacity * hytank.min_filling_ratio, hytank.capacity * hytank.max_filling_ratio, nloh)
    cycles_grid = range(0, bt.lifetime_cycles / 5, ncycles)
    RT_grid = range(0, 9000, nRT)



    #  dbt=0.05*bt.energy_rated
    #  dh2=0.05*fc.power_rated

    
    IBS = BSpline(Linear())

    nw = 11
    w_law = Normal(0, σ)
    demand = range(-3σ, 3σ, nw)
    prob_demand = pdf.(w_law, demand) / (sum(pdf(w_law, demand)))

    Jopt = zeros(Float32, nPnl, nloh, nEbt, ncycles, nRT, nRT, K + 1)
    # Policy over time
    Pbt_grid = zeros(nthreads(), 25)

    policy = zeros(Float16, nPnl, nloh, nEbt, ncycles, nRT, nRT, K)
    polh2 = zeros(Float16, nPnl, nloh, nEbt, ncycles, nRT, nRT, K)
    u_rb = zeros(Float16, nPnl, nloh, nEbt, ncycles, nRT, nRT, 2)

    @inbounds @threads for pnl = 1:nPnl

        @inbounds for loh = 1:nloh, ebt = 1:nEbt, cy = 1:ncycles, rt1 = 1:nRT, rt2 = 1:nRT
            xk = State(Pnl_grid[pnl], LoH_grid[loh], Ebt_grid[ebt], cycles_grid[cy], RT_grid[rt1], RT_grid[rt2])
            Jopt[pnl, loh, ebt, cy, rt1, rt2, K+1] = term(xk, K, mg) #arranger ça
            u_rb[pnl, loh, ebt, cy, rt1, rt2, :] .= dispatch_1(xk, mg)[:]
        end
    end


    @inbounds for k = K:-1:1

        Jopt_next = Interpolations.scale(interpolate!(selectdim(Jopt, ndims(Jopt), k + 1), IBS), Pnl_grid, LoH_grid, Ebt_grid, cycles_grid, RT_grid, RT_grid)

        @inbounds @threads for pnl = 1:nPnl

            @inbounds for loh = 1:nloh, ebt = 1:nEbt, cy = 1:ncycles, rt1 = 1:nRT, rt2 = 1:nRT

                xk = State(Pnl_grid[pnl], LoH_grid[loh], Ebt_grid[ebt], cycles_grid[cy], RT_grid[rt1], RT_grid[rt2])

                nu = u_range!(selectdim(Pbt_grid, 1, threadid()), xk, mg, 25, 12, 12) #peut changer de taille à régler
                err=(xk.Pnl-ref[k])
                #CREATE UGRID
                Jk_xu_min = Inf
                u_opt = Pbt_grid[threadid(), 1]
                ph_opt = 0.0

                @inbounds for u_i = 1:nu  #gérer les situations où il n'y pas de u
                    uk = Pbt_grid[threadid(), u_i]
                    Jk_xu = 0.0

                    PH2 = 0.0
                     x_next, Pshed, PH2 = dynamic(xk, uk, 0.0, mg)

                    @inbounds @fastmath for iw in eachindex(demand) # for each uncertainty value                
                        Jk_xu += Jopt_next(ref[k+1] + (err * ϕ + demand[iw]) , x_next.LoH, x_next.Ebatt, x_next.cycles, x_next.RT_fc, x_next.RT_el) * prob_demand[iw]

                    end # each uncertainty

                    Jk_xu += pen *max(0.0, Pshed)
                    if (Jk_xu < Jk_xu_min)
                        u_opt = uk
                        Jk_xu_min = Jk_xu
                        ph_opt = PH2
                    end

                end #each control


                Jopt[pnl, loh, ebt, cy, rt1, rt2, k] = Jk_xu_min
                policy[pnl, loh, ebt, cy, rt1, rt2, k] = u_opt
                polh2[pnl, loh, ebt, cy, rt1, rt2, k] = ph_opt

            end # each state
        end
        #      
    end #each instant


    return Jopt, policy, polh2, u_rb

end


function DPrecursion_jopt(K::Int64, pen::Float64, max_pnl::Float64, min_pnl::Float64,ref::Vector{Float64},σ::Float64,ϕ::Float64, mg)

    bt = mg.storage
    hytank = mg.tanks.h2Tank
    dt = mg.project.timestep
    nPnl = 16
    nEbt = 8
    nloh = 8
    nRT = 10
    ncycles = 10

    Pnl_grid = range(min_pnl, max_pnl, nPnl)
    Ebt_grid = range(bt.energy_rated * bt.SoC_min, bt.energy_rated * bt.SoC_max, nEbt)
    LoH_grid = range(hytank.capacity * hytank.min_filling_ratio, hytank.capacity * hytank.max_filling_ratio, nloh)
    cycles_grid = range(0, bt.lifetime_cycles / 5, ncycles)
    RT_grid = range(0, 9000, nRT)


    
    IBS = BSpline(Linear())

    nw = 11
    w_law = Normal(0, σ)
    demand = range(-3σ, 3σ, nw)
    prob_demand = pdf.(w_law, demand) / (sum(pdf(w_law, demand)))

    Jopt = zeros(Float32, nPnl, nloh, nEbt, ncycles, nRT, nRT, K + 1)
    # Policy over time
    Pbt_grid = zeros(nthreads(), 25)

    @inbounds @threads for pnl = 1:nPnl

        @inbounds for loh = 1:nloh, ebt = 1:nEbt, cy = 1:ncycles, rt1 = 1:nRT, rt2 = 1:nRT
            xk = State(Pnl_grid[pnl], LoH_grid[loh], Ebt_grid[ebt], cycles_grid[cy], RT_grid[rt1], RT_grid[rt2])
            Jopt[pnl, loh, ebt, cy, rt1, rt2, K+1] = term(xk, K, mg) #arranger ça
        end
    end


    @inbounds for k = K:-1:1

        Jopt_next = Interpolations.scale(interpolate!(selectdim(Jopt, ndims(Jopt), k + 1), IBS), Pnl_grid, LoH_grid, Ebt_grid, cycles_grid, RT_grid, RT_grid)

        @inbounds @threads for pnl = 1:nPnl

            @inbounds for loh = 1:nloh, ebt = 1:nEbt, cy = 1:ncycles, rt1 = 1:nRT, rt2 = 1:nRT

                xk = State(Pnl_grid[pnl], LoH_grid[loh], Ebt_grid[ebt], cycles_grid[cy], RT_grid[rt1], RT_grid[rt2])

                nu = u_range!(selectdim(Pbt_grid, 1, threadid()), xk, mg, 25, 12, 12) #peut changer de taille à régler
                err=(xk.Pnl-ref[k])
                #CREATE UGRID
                Jk_xu_min = Inf
                u_opt = Pbt_grid[threadid(), 1]


                @inbounds for u_i = 1:nu  #gérer les situations où il n'y pas de u
                    uk = Pbt_grid[threadid(), u_i]
                    Jk_xu = 0.0


                    x_next, Pshed, PH2 = dynamic(xk, uk, 0.0, mg)

                    @inbounds @fastmath for iw in eachindex(demand) # for each uncertainty value                
                        Jk_xu += Jopt_next(ref[k+1] + (err * ϕ + demand[iw]), x_next.LoH, x_next.Ebatt, x_next.cycles, x_next.RT_fc, x_next.RT_el) * prob_demand[iw]

                    end # each uncertainty

                    Jk_xu += pen * max(0.0, Pshed)
                    
                    if (Jk_xu < Jk_xu_min)
                        u_opt = uk
                        Jk_xu_min = Jk_xu

                    end

                end #each control


                Jopt[pnl, loh, ebt, cy, rt1, rt2, k] = Jk_xu_min
            end # each state
        end
        #      

    end #each instant


    return Jopt

end

function DPrecursion_jopt_mc(K::Int64, pen::Float64, time_init::Vector{Float64}, mi_l::Vector{Int64}, Pnl_cate, Pnl_prob, mg)

    bt = mg.storage
    # el=mg.electrolyzer[1]
    #fc= mg.dispatchables.fuel_cell[1]
    hytank = mg.tanks.h2Tank
    dt = mg.project.timestep
    nPnl = 16
    nEbt = 8
    nloh = 8
    nRT = 10
    ncycles = 10

    #Pnl_grid=range(min_pnl,max_pnl,nPnl)
    n_cat = length(Pnl_cate[1, 1, :]) - 1

    Pnl_grid = zeros(Float64, 12, 24, nPnl)
    for i = 1:12, j = 1:24
        Pnl_grid[i, j, :] = range(Pnl_cate[i, j, 1], Pnl_cate[i, j, end], nPnl)
    end

    Ebt_grid = range(bt.energy_rated * bt.SoC_min, bt.energy_rated * bt.SoC_max, nEbt)
    LoH_grid = range(hytank.capacity * hytank.min_filling_ratio, hytank.capacity * hytank.max_filling_ratio, nloh)
    cycles_grid = range(0, bt.lifetime_cycles / 5, ncycles)
    RT_grid = range(0, 9000, nRT)



    #  dbt=0.05*bt.energy_rated
    #  dh2=0.05*fc.power_rated

    time = time_init

    IBS = BSpline(Linear())



    Jopt = zeros(Float32, nPnl, nloh, nEbt, ncycles, nRT, nRT, K + 1)
    # Policy over time
    Pbt_grid = zeros(nthreads(), 25)

    @inbounds @threads for pnl = 1:nPnl

        @inbounds for loh = 1:nloh, ebt = 1:nEbt, cy = 1:ncycles, rt1 = 1:nRT, rt2 = 1:nRT
            xk = State(Pnl_grid[12, 24, pnl], LoH_grid[loh], Ebt_grid[ebt], cycles_grid[cy], RT_grid[rt1], RT_grid[rt2])
            Jopt[pnl, loh, ebt, cy, rt1, rt2, K+1] = term(xk, K, mg) #arranger ça
            
        end
    end


    @inbounds for k = K:-1:1
        time_next = time
        time_dyn_backward!(time, dt, mi_l)

        m = Int(time[3])
        h = Int(time[1])
        m_next = Int(time_next[3])
        h_next = Int(time_next[1])
        Jopt_next = Interpolations.scale(interpolate!(selectdim(Jopt, ndims(Jopt), k + 1), IBS), range(Pnl_cate[m_next, h_next+1, 1], Pnl_cate[m_next, h_next+1, end], nPnl), LoH_grid, Ebt_grid, cycles_grid, RT_grid, RT_grid)


        @inbounds  @threads for pnl = 1:nPnl

            @inbounds for loh = 1:nloh, ebt = 1:nEbt, cy = 1:ncycles, rt1 = 1:nRT, rt2 = 1:nRT

                xk = State(Pnl_grid[m, h+1, pnl], LoH_grid[loh], Ebt_grid[ebt], cycles_grid[cy], RT_grid[rt1], RT_grid[rt2])
                cat = Int(div(xk.Pnl - Pnl_cate[m, h+1, 1], Pnl_cate[m, h+1, 2] - Pnl_cate[m, h+1, 1]) + 1)
                if cat > 16
                     cat = 16
                end
                nu = u_range!(selectdim(Pbt_grid, 1, threadid()), xk, mg, 25, 12, 12) #peut changer de taille à régler

                #CREATE UGRID
                Jk_xu_min = Inf32
                u_opt = Pbt_grid[threadid(), 1]


                @inbounds for u_i = 1:nu  #gérer les situations où il n'y pas de u
                    uk = Pbt_grid[threadid(), u_i]
                    Jk_xu = 0.0
                    x_next, Pshed, PH2 = dynamic(xk, uk, 0.0, mg)

                   @inbounds @fastmath  for iw = 1:n_cat # for each uncertainty value  
                        if Pnl_prob[m, h+1, iw, cat] > 1e-10
                            Jk_xu += Jopt_next(Pnl_cate[m_next, h_next+1, iw+1] , x_next.LoH, x_next.Ebatt, x_next.cycles, x_next.RT_fc, x_next.RT_el) * Pnl_prob[m, h+1, iw, cat]
                            @assert Jk_xu >= 0. " problem, $m, $(h+1), $cat, $Jk_xu, $(Pnl_prob[m, h+1, iw, cat])"
                        end
                    end # each uncertainty

                    Jk_xu += pen * max(0.0, Pshed)
                    if (Jk_xu < Jk_xu_min)
                        u_opt = uk
                        Jk_xu_min = Jk_xu

                    end

                end #each control

                @assert Jk_xu_min >= Jopt[pnl, loh, ebt, cy, rt1, rt2, K+1] - 0.1 " problem, $m, $(h+1), $cat,$xk, $(dynamic(xk, u_opt, 0.0, mg)), $u_opt,$( Jk_xu_min), $(Jopt[pnl,loh,ebt,cy,rt1,rt2,K+1])"

                Jopt[pnl, loh, ebt, cy, rt1, rt2, k] = Jk_xu_min


            end # each state
        end
        #      


    end #each instant


    return Jopt

end

function DPrecursion_uopt(K::Int64, pen::Float64,  max_pnl::Float64, min_pnl::Float64, ref,σ,ϕ, mg)

    bt = mg.storage
    # el=mg.electrolyzer[1]
    #fc= mg.dispatchables.fuel_cell[1]
    hytank = mg.tanks.h2Tank
    dt = mg.project.timestep
    nPnl = 16
    nEbt = 8
    nloh = 8
    nRT = 10
    ncycles = 10

    Pnl_grid = range(min_pnl, max_pnl, nPnl)
    Ebt_grid = range(bt.energy_rated * bt.SoC_min, bt.energy_rated * bt.SoC_max, nEbt)
    LoH_grid = range(hytank.capacity * hytank.min_filling_ratio, hytank.capacity * hytank.max_filling_ratio, nloh)
    cycles_grid = range(0, bt.lifetime_cycles / 5, ncycles)
    RT_grid = range(0, 9000, nRT)



    #  dbt=0.05*bt.energy_rated
    #  dh2=0.05*fc.power_rated

    
    IBS = BSpline(Linear())

    nw = 11
    w_law = Normal(0, σ)
    demand = range(-3σ, 3σ, nw)
    prob_demand = pdf.(w_law, demand) / (sum(pdf(w_law, demand)))

    Jopt = zeros(Float32, nPnl, nloh, nEbt, ncycles, nRT, nRT, K + 1)
    # Policy over time
    Pbt_grid = zeros(nthreads(), 25)

    policy = zeros(Float16, nPnl, nloh, nEbt, ncycles, nRT, nRT, K)

    u_rb = zeros(Float16, nPnl, nloh, nEbt, ncycles, nRT, nRT, 2)

    @inbounds @threads for pnl = 1:nPnl

        @inbounds for loh = 1:nloh, ebt = 1:nEbt, cy = 1:ncycles, rt1 = 1:nRT, rt2 = 1:nRT
            xk = State(Pnl_grid[pnl], LoH_grid[loh], Ebt_grid[ebt], cycles_grid[cy], RT_grid[rt1], RT_grid[rt2])
            Jopt[pnl, loh, ebt, cy, rt1, rt2, K+1] = term(xk, K, mg) #arranger ça
            u_rb[pnl, loh, ebt, cy, rt1, rt2, :] .= dispatch_1(xk, mg)[:]
        end
    end


    @inbounds for k = K:-1:1

        Jopt_next = Interpolations.scale(interpolate!(selectdim(Jopt, ndims(Jopt), k + 1), IBS), Pnl_grid, LoH_grid, Ebt_grid, cycles_grid, RT_grid, RT_grid)

        @inbounds @threads for pnl = 1:nPnl

            @inbounds for loh = 1:nloh, ebt = 1:nEbt, cy = 1:ncycles, rt1 = 1:nRT, rt2 = 1:nRT

                xk = State(Pnl_grid[pnl], LoH_grid[loh], Ebt_grid[ebt], cycles_grid[cy], RT_grid[rt1], RT_grid[rt2])
                 err=(xk.Pnl-ref[k])
                nu = u_range!(selectdim(Pbt_grid, 1, threadid()), xk, mg, 25, 12, 12) #peut changer de taille à régler

                #CREATE UGRID
                Jk_xu_min = Inf
                u_opt = Pbt_grid[threadid(), 1]


                @inbounds for u_i = 1:nu  #gérer les situations où il n'y pas de u
                    uk = Pbt_grid[threadid(), u_i]
                    Jk_xu = 0.0


                     x_next, Pshed, PH2 = dynamic(xk, uk, 0.0, mg)

                    @inbounds @fastmath for iw in eachindex(demand) # for each uncertainty value                
                        Jk_xu += Jopt_next(ref[k+1] + (err * ϕ + demand[iw]) , x_next.LoH, x_next.Ebatt, x_next.cycles, x_next.RT_fc, x_next.RT_el) * prob_demand[iw]

                    end # each uncertainty

                    Jk_xu += pen * max(0.0, Pshed)
                    if (Jk_xu < Jk_xu_min)
                        u_opt = uk
                        Jk_xu_min = Jk_xu

                    end

                end #each control


                Jopt[pnl, loh, ebt, cy, rt1, rt2, k] = Jk_xu_min
                policy[pnl, loh, ebt, cy, rt1, rt2, k] = u_opt

            end # each state
        end
        #      

       

    end #each instant


    return policy, u_rb

end


function DPrecursion_rb_cost(K::Int64, pen::Float64, max_pnl::Float64, min_pnl::Float64,ref,σ,ϕ, mg)

    bt = mg.storage
   
    hytank = mg.tanks.h2Tank
    dt = mg.project.timestep
    nPnl = 16
    nEbt = 8
    nloh = 8
    nRT = 10
    ncycles = 10

    Pnl_grid = range(min_pnl, max_pnl, nPnl)
    Ebt_grid = range(bt.energy_rated * bt.SoC_min, bt.energy_rated * bt.SoC_max, nEbt)
    LoH_grid = range(hytank.capacity * hytank.min_filling_ratio, hytank.capacity * hytank.max_filling_ratio, nloh)
    cycles_grid = range(0, bt.lifetime_cycles / 5, ncycles)
    RT_grid = range(0, 9000, nRT)


    IBS = BSpline(Linear())

    nw = 11
     w_law = Normal(0, σ)
    demand = range(-3σ, 3σ, nw)
    prob_demand = pdf.(w_law, demand) / (sum(pdf(w_law, demand)))

    Jopt = zeros(Float32, nPnl, nloh, nEbt, ncycles, nRT, nRT, K + 1)
    # Policy over time




    u_rb = zeros(Float16, nPnl, nloh, nEbt, ncycles, nRT, nRT, 2)

    @inbounds @threads for pnl = 1:nPnl

        @inbounds for loh = 1:nloh, ebt = 1:nEbt, cy = 1:ncycles, rt1 = 1:nRT, rt2 = 1:nRT
            xk = State(Pnl_grid[pnl], LoH_grid[loh], Ebt_grid[ebt], cycles_grid[cy], RT_grid[rt1], RT_grid[rt2])
            Jopt[pnl, loh, ebt, cy, rt1, rt2, K+1] = term(xk, K, mg) #arranger ça
            u_rb[pnl, loh, ebt, cy, rt1, rt2, :] .= dispatch_1(xk, mg)[:]
        end
    end


    @inbounds for k = K:-1:1

        Jopt_next = Interpolations.scale(interpolate!(selectdim(Jopt, ndims(Jopt), k + 1), IBS), Pnl_grid, LoH_grid, Ebt_grid, cycles_grid, RT_grid, RT_grid)

        @inbounds @threads for pnl = 1:nPnl

            @inbounds for loh = 1:nloh, ebt = 1:nEbt, cy = 1:ncycles, rt1 = 1:nRT, rt2 = 1:nRT

                xk = State(Pnl_grid[pnl], LoH_grid[loh], Ebt_grid[ebt], cycles_grid[cy], RT_grid[rt1], RT_grid[rt2])
                 err=(xk.Pnl-ref[k])




                Jk_xu = 0.0
                uk = u_rb[pnl, loh, ebt, cy, rt1, rt2, 1]

                x_next, Pshed, PH2 = dynamic(xk, uk, 0.0, mg)

                @inbounds @fastmath for iw in eachindex(demand) # for each uncertainty value                
                Jk_xu += Jopt_next(ref[k+1] + (err * ϕ + demand[iw]) , x_next.LoH, x_next.Ebatt, x_next.cycles, x_next.RT_fc, x_next.RT_el) * prob_demand[iw]
                end # each uncertainty

                Jk_xu += pen * max(0.0, Pshed)



                Jopt[pnl, loh, ebt, cy, rt1, rt2, k] = Jk_xu


            end # each state
        end
        #      

    end #each instant


    return Jopt[:, :, :, :, :, :, 1]

end
function sdp_simu_jopt(Jopt, n, pen,ϕ,σ, mg1,ref::Vector{Float64},max_pnl::Float64, min_pnl::Float64,Pnl::Vector{Float64})
    Pren = zeros(n)
    Pren = production(mg1.nondispatchables[1]) .+ production(mg1.nondispatchables[2])
    Pbt = zeros(n)
    Pel = zeros(n)
    Pfc = zeros(n)
    Ebt = zeros(n + 1)
    cycles = 0.0
    Loh = zeros(n + 1)
    Lof = zeros(n + 1)
    Pgen = zeros(n)
    Phb = zeros(n)
    Pshed = zeros(n)
    Pdump = zeros(n)
    Pspill = zeros(n)
   
    IBS = BSpline(Linear())
   
    Ebt[1] = 0.
    Loh[1] = mg1.tanks.h2Tank.ini_filling_ratio * mg1.tanks.h2Tank.capacity
    RT_fc, RT_el = 0, 0
    nEbt = size(Jopt)[3]
    nPnl = size(Jopt)[1]
    nloh = size(Jopt)[2]
    nRT = size(Jopt)[5]
    ncycles = size(Jopt)[4]
    Ebt_grid = range(bt.energy_rated * bt.SoC_min, bt.energy_rated * bt.SoC_max, nEbt)
    Pnl_grid = range(min_pnl, max_pnl, nPnl)
    LoH_grid = range(hytank.capacity * hytank.min_filling_ratio, hytank.capacity * hytank.max_filling_ratio, nloh)
    RT_grid = range(0, 9000, nRT)
    cycles_grid = range(0, bt.lifetime_cycles / 5, ncycles)
    Pbt_grid = zeros(50)
    nw = 11
    w_law = Normal(0, σ)
    demand = range(-3σ, 3σ, nw)
    prob_demand = pdf.(w_law, demand) / (sum(pdf(w_law, demand)))
    for k = 1:n
        J_interp = Interpolations.extrapolate(Interpolations.scale(Interpolations.interpolate!(selectdim(Jopt, ndims(Jopt), k + 1), IBS), Pnl_grid, LoH_grid, Ebt_grid, cycles_grid, RT_grid, RT_grid), Line())

        xk = State(Pnl[k], Loh[k], Ebt[k], cycles, RT_fc, RT_el)
        err=(xk.Pnl-ref[k])
        nu = u_range!(Pbt_grid, xk, mg1, 20, 20, 10)
       
        
        Jk_xu_min = Inf
        u_opt = 0.0
        ph_opt = 0.0

        @inbounds for u_i = 1:nu  #gérer les situations où il n'y pas de u
            uk = Pbt_grid[u_i]
            Jk_xu = 0.0

           x_next, Ps, Ph2 = dynamic(xk, uk, 0.0, mg1)

                @inbounds @fastmath for iw in eachindex(demand) # for each uncertainty value                
                Jk_xu += J_interp(ref[k+1] + (err * ϕ + demand[iw]) , x_next.LoH, x_next.Ebatt, x_next.cycles, x_next.RT_fc, x_next.RT_el) * prob_demand[iw]
                end # each uncertainty



            Jk_xu += pen * max(0.0, Ps)
            if (Jk_xu < Jk_xu_min)
                u_opt = uk
                Jk_xu_min = Jk_xu
                ph_opt = Ph2
            end

        end #each control
        x_next, Ps, Ph2 = dynamic(xk, u_opt, 0.0, mg1)
        Pbt[k] = u_opt
        Pel[k] = ifelse(Ph2 > 0.0, 0.0, -Ph2)
        Pfc[k] = ifelse(Ph2 > 0.0, Ph2, 0.0)
        Ebt[k+1] = x_next.Ebatt
        Loh[k+1] = x_next.LoH
        cycles = x_next.cycles
        Pspill[k] = ifelse(Ps > 0.0, 0.0, -Ps)
        Pshed[k] = ifelse(Ps > 0.0, Ps, 0.0)
        RT_el=x_next.RT_el
        RT_fc=x_next.RT_fc

    end
    traj = OperationTraj(Pnl, Pshed, Pren, Pgen, Pfc, Pel, Phb, Ebt, Pbt, Loh, Lof, Pspill, Pdump)
    stats = aggregation(mg1, traj)
    # Eval the microgrid costs
    costs = economics(mg1, stats)

    return traj, stats, costs

end


function sdp_simu_jopt_mc(Jopt, n, pen, mg1, Pnl_cate, Pnl_prob,Pnl)
    Pren = zeros(n)
   Pren = production(mg1.nondispatchables[1]) .+ production(mg1.nondispatchables[2])

    Pbt = zeros(n)
    Pel = zeros(n)
    Pfc = zeros(n)
    Ebt = zeros(n + 1)
    cycles = 0.0
    Loh = zeros(n + 1)
    Lof = zeros(n + 1)
    Pgen = zeros(n)
    Phb = zeros(n)
    Pshed = zeros(n)
    Pdump = zeros(n)
    Pspill = zeros(n)
    
    IBS = BSpline(Linear())
    
    Ebt[1] = 0.
    Loh[1] = mg1.tanks.h2Tank.ini_filling_ratio * mg1.tanks.h2Tank.capacity
    RT_fc, RT_el = 0, 0
    nEbt = size(Jopt)[3]
    nPnl = size(Jopt)[1]
    nloh = size(Jopt)[2]
    nRT = size(Jopt)[5]
    ncycles = size(Jopt)[4]
    Ebt_grid = range(bt.energy_rated * bt.SoC_min, bt.energy_rated * bt.SoC_max, nEbt)
    Pnl_grid = zeros(Float64, 12, 24, nPnl)
    n_cat = length(Pnl_cate[1, 1, :]) - 1


    for i = 1:12, j = 1:24
        Pnl_grid[i, j, :] = range(Pnl_cate[i, j, 1], Pnl_cate[i, j, end], nPnl)
    end

    LoH_grid = range(hytank.capacity * hytank.min_filling_ratio, hytank.capacity * hytank.max_filling_ratio, nloh)
    RT_grid = range(0, 9000, nRT)
    cycles_grid = range(0, bt.lifetime_cycles / 5, ncycles)
    Pbt_grid = zeros(100)

    for k = 1:n
        if k < n
            time_next = [Float64(data.Hour[k+1]), Float64(data.Day[k+1]), Float64(data.Month[k+1])]

        else
            time_next = [0., 1, 1]
        end

        time = [Float64(data.Hour[k]), Float64(data.Day[k]), Float64(data.Month[k])]
        m = Int(time[3])
        h = Int(time[1])
        m_next = Int(time_next[3])
        h_next = Int(time_next[1])


        J_interp = Interpolations.extrapolate(Interpolations.scale(Interpolations.interpolate!(selectdim(Jopt, ndims(Jopt), k + 1), IBS), range(Pnl_cate[m_next, h_next+1, 1], Pnl_cate[m_next, h_next+1, end], nPnl),
                LoH_grid, Ebt_grid, cycles_grid, RT_grid, RT_grid), Line())

        xk = State(Pnl[k], Loh[k], Ebt[k], cycles, RT_fc, RT_el)
        cat = Int(div(xk.Pnl - Pnl_cate[m, h+1, 1], Pnl_cate[m, h+1, 2] - Pnl_cate[m, h+1, 1]) + 1)

        if cat > 16
            cat = 16
        end
        
        if cat < 1
            cat = 1
        end
        
        @assert cat >= 1 && cat <= 16 "error $cat, $(xk.Pnl),  $(h+1), $m, $(Pnl_cate[m,h+1,1]), $(xk.Pnl-Pnl_cate[m,h+1,1]), $(Pnl_cate[m,h+1,2]-Pnl_cate[m,h+1,1]))"
        nu = u_range!(Pbt_grid, xk, mg1, 20, 50, 50)
        #  nu_bt= Int(ceil((umax.P_bt-umin.P_bt)/dbt)) + 1
        #  nu_h2= Int(ceil((umax.P_H2-umin.P_H2)/dh2)) + 1




        #CREATE UGRID
        Jk_xu_min = Inf
        u_opt = 0.0
        ph_opt = 0.0


        @inbounds for u_i = 1:nu  #gérer les situations où il n'y pas de u
            uk = Pbt_grid[u_i]
            Jk_xu = 0.0
            x_next, Ps, Ph2 = dynamic(xk, uk, 0.0, mg1)

            @fastmath for iw = 1:n_cat # for each uncertainty value  
                if Pnl_prob[m, h+1, iw, cat] != 0.0

                    Jk_xu += J_interp(Pnl_cate[m_next, h_next+1, iw+1], x_next.LoH, x_next.Ebatt, x_next.cycles, x_next.RT_fc, x_next.RT_el) * Pnl_prob[m, h+1, iw, cat]

                end
            end # each uncertainty

            #   @assert Jk_xu=0.0!=0 " cout nul"
            Jk_xu += pen * max(0.0, Ps)
            if (Jk_xu < Jk_xu_min)
                u_opt = uk
                Jk_xu_min = Jk_xu
                ph_opt = Ph2
            end

        end #each control
        x_next, Ps, Ph2 = dynamic(xk, u_opt, 0.0, mg1)
        Pbt[k] = u_opt
        Pel[k] = ifelse(Ph2 > 0.0, 0.0, -Ph2)
        Pfc[k] = ifelse(Ph2 > 0.0, Ph2, 0.0)
        Ebt[k+1] = x_next.Ebatt
        Loh[k+1] = x_next.LoH
        cycles = x_next.cycles
        RT_el=x_next.RT_el
        RT_fc=x_next.RT_fc
        Pspill[k] = ifelse(Ps > 0.0, 0.0, -Ps)
        Pshed[k] = ifelse(Ps > 0.0, Ps, 0.0)


    end
    traj = OperationTraj(Pnl, Pshed, Pren, Pgen, Pfc, Pel, Phb, Ebt, Pbt, Loh, Lof, Pspill, Pdump)
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
    m=Int64((n-1)/2)
    new_sig=zeros(length(signal)+2m)
    new_sig[1:m]=signal[end-m+1:end]
    new_sig[m+1:end-m]=signal[:]
    new_sig[end-m+1:end]=signal[1:m]

    return hma(new_sig,n)[m+1:end-m]

end
function plot_periodigram(signal)
    n=length(signal)
    w=range(1,n/2,n)./n
    i=zeros(Float64,n)
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
    nj=Int64(n/d)
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