# "Typical" parameters for a Microgrid project which contains
# - solar PV
# - Battery
# - generator
# and parameters for the Project itself
#
# Remark: time series are not included (Pload, Solar irradiance IT)

# Project
lifetime = 25 # yr
discount_rate = 0.05
timestep = 1. # h

# Diesel generator
power_rated_DG = 1800.
min_load_ratio = 0
F0 = 0.0
F1 = 0.240
fuel_cost = 1.
investiment_cost_DG = 400.
om_cost_DG = 0.02
replacement_cost_DG = 400.
salvage_cost_DG = 400.
lifetime_DG = 15000

# Battery energy storage
energy_initial = 0.
energy_max = 9000
energy_min = 0
power_min = -1.0*energy_max
power_max = +1.0*energy_max
loss = 0.05
investiment_cost_BT = 350.
om_cost_BT = 10.
replacement_cost_BT = 350.
salvage_cost_BT = 350.
lifetime_BT = 15
lifetime_thrpt = 3000

# Photovoltaic generation
power_rated_PV = 6000 # kW
fPV = 1.
IS = 1.0
investiment_cost_PV = 1200.
om_cost_PV = 20.
replacement_cost_PV = 1200.
salvage_cost_PV = 1200.
lifetime_PV = 25