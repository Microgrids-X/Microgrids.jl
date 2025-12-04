# Parameters and time series for a Microgrid project with *wind* and *solar* sources,
# plus a *battery* and a *dispatchable generator*.
#
# Values gathered from the Microgrid_Wind-Solar.ipynb notebook,
# to be used in the sizing optimization notebook.

using Microgrids # for capacity_from_wind
using CSV, DataFrames

println("Data definition for Microgrid with wind, solar and Hybrid Storage System...")

### Input time series
data = DataFrame(CSV.File("$(@__DIR__)/Ouessant_data_2016.csv"))

# Simulation steps
nsteps = length(data.Load)
@assert(nsteps == 8760) # 1 year at an hourly timestep

# Split load, solar and wind data:
Pload = Float32.(data.Load)# kW
Ppv1k =  Float32.(data.Ppv1k./ 1000); # convert to kW/kWp
wind_speed = data.Wind # m/s

# Calibrate wind speed data against a mast measurement
ws_gain = 1.059f0 # ratio of Mast's mean /PVGIS' mean
wind_speed = Float32.(ws_gain*wind_speed)

# Generate wind power capacity factor time series from wind speed
S_D52 =Float32(pi * (52/2)^2 ) # rotor swept area m²
TSP_D52 = Float32(900e3/S_D52) # W/m²
v_out = 25.0f0 # m/s
Cp_D52, α_D52 = 0.521f0, 3.1f0; # fitted from actual power curve
cf_wind = capacity_from_wind.(wind_speed; TSP=TSP_D52, Cp=Cp_D52, v_out=v_out, α=α_D52)

### Microgrid parameters

# Project parameters
lifetime = 25f0 # yr
discount_rate = 0.05f0
timestep = 1.f0 # h

# Parameters chosen to be common to all Components (but they could differ)
replacement_price_ratio = 1.0f0
salvage_price_ratio = 1.0f0

# Dispatchable generator (Diesel)
power_rated_gen = 2000.f0  # rated power (kW) → /2 to see some load shedding
# fuel parameters:
fuel_intercept = 0.0f0 # fuel curve intercept (L/h/kW_max)
fuel_slope = 0.240f0 # fuel curve slope (L/h/kW)
# prices:
fuel_price = 1.0f0 # fuel price ($/L)
investment_price_gen = 400.f0 # initial investiment price ($/kW)
om_price_hour_gen = 0.02f0 # operation & maintenance price ($/kW/h of operation)
om_price_gen = 0.f0 #
# other technical parameters
lifetime_gen_y = 25.f0
lifetime_gen_h = 20000.f0 # generator lifetime (h)
lifetime_gen_starts = 5000.f0 # generator lmax starts on
load_ratio_min_gen = 0.1f0 # minimum load ratio ∈ [0,1]
input_unit_gen = "L" # fuel quantity unit (used in fuel price and consumption curve parameters)
output_unit_gen = "KWh"

#fuel tank

capacity_rated_ftank = 0.f0 # rated power capacity (L)
investment_price_ftank = 0.0f0 # initial investment price  ($/L)
om_price_ftank = 0.f0# operation and maintenance price ($/L/y)
lifetime_ftank = 25.f0 # calendar lifetime (y)
loss_factor_ftank = 0.f0 # hydrogen used on site 
fuel_ini_ratio = 0.f0 # Initial load ratio ∈ [0,1]
fuel_min_ratio = 0.f0 # minimum load ratio ∈ [0,1]
fuel_max_ratio = 1.f0 # maximum load ratio ∈ [0,1]

#H2 Tank
capacity_rated_hytank = 12000.f0 # rated power capacity (kg)
investment_price_hytank = 400.f0 # initial investment price  ($/kg)
hy_price = 14.f0 # initial hydrogen price ($/kg)
om_price_hytank = 4.f0 # operation and maintenance price ($/kg/y)
lifetime_hytank = 25.f0 # calendar lifetime (y)
loss_factor_hytank = 0.f0 # hydrogen used on site 
LoH_ini_ratio = 0.0f0 # Initial load ratio ∈ [0,1]
LoH_min_ratio = 0.0f0 # minimum load ratio ∈ [0,1]
LoH_max_ratio = 1.f0 # maximum load ratio ∈ [0,1]


#Fuel cell
power_rated_fc = 1800.f0 # rated power capacity (KW)
cons_intercept_fc=0.f0 #
cons_rate_fc = 0.0625f0 # consumption rate (KgH2/KWhé)
cons_price_fc = 0.f0 #
investment_price_fc = 1600.f0 # initial investment price  ($/KW)
om_price_fc = 45.f0 # operation and maintenance price ($/kW/y)
om_price_hour_fc = 0.f0 # operation and maintenance price ($/kW/h)
lifetime_fc_y = 15.f0
lifetime_fc_h = 45000.f0 # Fuel Cell lifetime (h)
lifetime_fc_starts = 5000.f0# fuel cell maximum starts on
load_min_ratio_fc = 0.05f0 # minimum load ratio ∈ [0,1]
input_unit_fc= "Kg"
output_unit_fc="KWh"

#Haber Bosch numbers to update!!! they have been taken arbitrary
power_rated_hb = 4000.f0 # rated power capacity (KW)
cons_intercept_hb=0.f0 #
cons_rate_hb = 0.344f0 # consumption rate (KWhé/KgH2)
cons_price_hb = 0.f0 #
investment_price_hb_only = Float32(0.418*580) #energy consumption of HB process for 1 kW of HB+ASU multiplied by CAPEX of HB ($/kW) [Achour et al., 2025]
investment_price_asu_only =Float32( 0.582*224) #energy consumption of ASU for 1 kW of HB+ASU multiplied by CAPEX of ASU ($/kW) [Achour et al., 2025]
investment_price_hb = investment_price_asu_only+investment_price_hb_only # initial investment price ($/KW)
om_price_hb =Float32(0.017*1000) # operation and maintenance price ($/kW/y)[Achour et al., 2025]
om_price_hour_hb = 0.f0 # operation and maintenance price ($/kW/h)
lifetime_hb_y = 25.f0
lifetime_hb_h = 45000.f0 # FOR NOW NOT USED. THIS IS THE VALUE FOR FUELL CELL.
lifetime_hb_starts = 5000.f0 # fuel cell maximum starts on
load_min_ratio_hb= 0.3f0 # minimum load ratio ∈ [0,1]
input_unit_hb= "Kwh"
output_unit_hb="Kg"



# Battery energy storage
energy_rated_sto = 5000.f0 # rated energy capacity (kWh)
investment_price_sto = 350.f0 # initial investiment price ($/kWh)
om_price_sto =10.f0 # operation and maintenance price ($/kWh/y)
lifetime_sto = 15.f0 # calendar lifetime (y)
lifetime_cycles = 3000.f0 # maximum number of cycles over life (1)
# Parameters with default values
charge_rate = 1.0f0 # max charge power for 1 kWh (kW/kWh = h^-1)
discharge_rate = 1.0f0 # max discharge power for 1 kWh (kW/kWh = h^-1)
loss_factor_sto = 0.05f0 # linear loss factor α (round-trip efficiency is about 1 − 2α) ∈ [0,1]
SoC_min = 0.f0 # minimum State of Charge ∈ [0,1]
SoC_ini = 0.f0 # initial State of Charge ∈ [0,1]
SoC_max = 1.f0 # initial State of Charge ∈ [0,1]

#electrolyzer
power_rated_elyz = 2000.f0 # rated power capacity (kW)
cons_intercept_elyz= 0.f0 # consumption rate (KWhé/KgH2)
cons_rate_elyz = 56.f0  # consumption rate (KWhé/KgH2)
cons_price_elyz = 0.f0#
investment_price_elyz = 1600.f0 # initial investment price  ($/kW)
om_price_elyz = 44.f0 # operation and maintenance price ($/kW/y)
om_price_hour_elyz = 0.f0 # operation and maintenance price ($/kW/h)
lifetime_elyz_y = 20.f0
lifetime_elyz_h = 45000.f0 #Electrolyzer lifetime (h)
lifetime_elyz_starts = 5000.f0 #Electrolyzer maximum starts on 
load_min_ratio_elyz = 0.05f0  # minimum load ratio ∈ [0,1]

input_unit_elyz= "KWh"
output_unit_elyz="Kg"

# Photovoltaic (PV) generation
power_rated_pv = 5000.f0 # rated power (kW)
irradiance = Ppv1k # global solar irradiance incident on the PV array (kW/m²)
investment_price_pv = 1200.f0 # initial investiment price ($/kW)
om_price_pv = 20.f0# operation and maintenance price ($/kW/y)
lifetime_pv = 25.f0 # lifetime (y)
# Parameters with default values
derating_factor_pv = 1.0f0 # derating factor (or performance ratio) ∈ [0,1]

# Wind power generation
power_rated_wind = 1800.f0 # rated power (kW)
investment_price_wind = 3500.f0 # initial investiment price ($/kW)
om_price_wind = 100.f0# operation and maintenance price ($/kW/y)
lifetime_wind = 25.f0 # lifetime (y)

default_sizing= Sizing( 0.0f0, 5000.f0,  5000.f0, 1800.f0,  2000.f0, 1800.f0,  12000.f0, 0.0f0,0.0f0)
capex_def=[400.f0, 0.0f0, 1600.f0, 500.f0, 1600.f0, 350.f0, 1200.f0, 3500.f0]
ini_filling_state=[0.0f0,0.0f0,0.0f0];
Eol=0.9
deg_ratio_rt=Float32(1.5e-6/0.7)# 10 μV*/h same ratio as Morin
deg_ratio_st=Float32(14e-6/0.7) # 106 μV*/h same ratio as Morin