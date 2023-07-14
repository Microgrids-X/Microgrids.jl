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

# Parameters common to all Components
replacement_price_ratio = 1.0
salvage_price_ratio = 1.0

# Dispatchable generator
power_rated_gen = 1800.
fuel_intercept = 0.0
fuel_slope = 0.240
fuel_price = 1.0
investment_price_gen = 400.
om_price_gen = 0.02
lifetime_gen = 15000.
load_ratio_min = 0.0
fuel_unit = "L"

# Battery energy storage
energy_rated = 9000.
investment_price_sto = 350.
om_price_sto = 10.
lifetime_sto = 15.
lifetime_cycles = 3000.
charge_rate = 1.0
discharge_rate = 1.0
loss_factor_sto = 0.05
SoC_min = 0.
SoC_ini = 0.

# Photovoltaic generation
power_rated_pv = 6000.
investment_price_pv = 1200.
om_price_pv = 20.
lifetime_pv = 25.
derating_factor_pv = 1.0