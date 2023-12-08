# Parameters and time series for a Microgrid project with *wind* and *solar* sources,
# plus a *battery* and a *dispatchable generator*.
#
# Values gathered from the Microgrid_Wind-Solar.ipynb notebook,
# to be used in the sizing optimization notebook.

using Microgrids # for capacity_from_wind
using CSV, DataFrames

println("Data definition for Microgrid with wind, solar, storage and generator...")

### Input time series
data = DataFrame(CSV.File("$(@__DIR__)/data/Ouessant_data_2016.csv"))

# Simulation steps
nsteps = length(data.Load)
@assert(nsteps == 8760) # 1 year at an hourly timestep

# Split load, solar and wind data:
Pload = data.Load # kW
Ppv1k =  data.Ppv1k ./ 1000; # convert to kW/kWp
wind_speed = data.Wind; # m/s

# Calibrate wind speed data against a mast measurement
ws_gain = 1.059 # ratio of Mast's mean /PVGIS' mean
wind_speed = ws_gain*wind_speed

# Generate wind power capacity factor time series from wind speed
S_D52 = pi * (52/2)^2 # rotor swept area m²
TSP_D52 = 900e3/S_D52 # W/m²
v_out = 25.0 # m/s
Cp_D52, α_D52 = 0.521, 3.1; # fitted from actual power curve
cf_wind = capacity_from_wind.(wind_speed; TSP=TSP_D52, Cp=Cp_D52, v_out=v_out, α=α_D52)

### Microgrid parameters

# Project parameters
lifetime = 25 # yr
discount_rate = 0.05
timestep = 1. # h

# Parameters chosen to be common to all Components (but they could differ)
replacement_price_ratio = 1.0
salvage_price_ratio = 1.0

# Dispatchable generator (Diesel)
power_rated_gen = 2000.  # rated power (kW) → /2 to see some load shedding
# fuel parameters:
fuel_intercept = 0.0 # fuel curve intercept (L/h/kW_max)
fuel_slope = 0.240 # fuel curve slope (L/h/kW)
# prices:
fuel_price = 1. # fuel price ($/L)
investment_price_gen = 400. # initial investiment price ($/kW)
om_price_gen = 0.02 # operation & maintenance price ($/kW/h of operation)
# other technical parameters
lifetime_gen = 15000. # generator lifetime (h)
load_ratio_min_gen = 0.0 # minimum load ratio ∈ [0,1]
input_unit_gen = "L" # fuel quantity unit (used in fuel price and consumption curve parameters)
output_unit_gen = "KWh"

#fuel tank

capacity_rated_ftank = 0. # rated power capacity (L)
investment_price_ftank = 0.0 # initial investment price  ($/L)
om_price_ftank = 0.# operation and maintenance price ($/L/y)
lifetime_ftank = 25. # calendar lifetime (y)
loss_factor_ftank = 0. # hydrogen used on site 
fuel_ini_ratio = 0. # Initial load ratio ∈ [0,1]
fuel_min_ratio = 0. # minimum load ratio ∈ [0,1]
fuel_max_ratio = 1. # maximum load ratio ∈ [0,1]

#H2 Tank
capacity_rated_hytank = 10000. # rated power capacity (kg)
investment_price_hytank = 0.8 # initial investment price  ($/kg)
om_price_hytank = 0.0025 # operation and maintenance price ($/kg/y)
lifetime_hytank = 25. # calendar lifetime (y)
loss_factor_hytank = 0. # hydrogen used on site 
LoH_ini_ratio = 0.2 # Initial load ratio ∈ [0,1]
LoH_min_ratio = 0.1 # minimum load ratio ∈ [0,1]
LoH_max_ratio = 1. # maximum load ratio ∈ [0,1]


#Fuel cell
power_rated_fc = 1800. # rated power capacity (KW)
cons_intercept_fc=0. #
cons_rate_fc = 0.0625 # consumption rate (KgH2/KWhé)
cons_price_fc = 0. #
investment_price_fc = 1000. # initial investment price  ($/KW)
om_price_fc = 10. # operation and maintenance price ($/kW/y)
lifetime_fc = 45000. # Fuel Cell lifetime (h)
load_min_ratio_fc = 0.05 # minimum load ratio ∈ [0,1]
input_unit_fc= "Kg"
output_unit_fc="KWh"


# Battery energy storage
energy_rated_sto = 5000. # rated energy capacity (kWh)
investment_price_sto = 350. # initial investiment price ($/kWh)
om_price_sto = 10. # operation and maintenance price ($/kWh/y)
lifetime_sto = 15. # calendar lifetime (y)
lifetime_cycles = 3000. # maximum number of cycles over life (1)
# Parameters with default values
charge_rate = 1.0 # max charge power for 1 kWh (kW/kWh = h^-1)
discharge_rate = 1.0 # max discharge power for 1 kWh (kW/kWh = h^-1)
loss_factor_sto = 0.05 # linear loss factor α (round-trip efficiency is about 1 − 2α) ∈ [0,1]
SoC_min = 0. # minimum State of Charge ∈ [0,1]
SoC_ini = 0. # initial State of Charge ∈ [0,1]
SoC_max = 1. # initial State of Charge ∈ [0,1]

#electrolyzer
power_rated_elyz = 2000. # rated power capacity (Kw)
cons_intercept_elyz= 0. # consumption rate (KWhé/KgH2)
cons_slope_elyz = 56.  # consumption rate (KWhé/KgH2)
cons_price_elyz = 0.#
investment_price_elyz = 1000. # initial investment price  ($/kW)
om_price_elyz = 10. # operation and maintenance price ($/kW/y)
lifetime_elyz = 45000. #Electrolyzer lifetime (h)
load_min_ratio_elyz = 0.05  # minimum load ratio ∈ [0,1]

input_unit_elyz= "KWh"
output_unit_elyz="Kg"

# Photovoltaic (PV) generation
power_rated_pv = 3000. # rated power (kW)
irradiance = Ppv1k # global solar irradiance incident on the PV array (kW/m²)
investment_price_pv = 1200. # initial investiment price ($/kW)
om_price_pv = 20.# operation and maintenance price ($/kW/y)
lifetime_pv = 25. # lifetime (y)
# Parameters with default values
derating_factor_pv = 1.0 # derating factor (or performance ratio) ∈ [0,1]

# Wind power generation
power_rated_wind = 900. # rated power (kW)
investment_price_wind = 3500. # initial investiment price ($/kW)
om_price_wind = 100.# operation and maintenance price ($/kW/y)
lifetime_wind = 25. # lifetime (y)

