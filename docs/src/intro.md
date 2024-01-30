# Microgrids.jl
 ### Tests the objective functions

Sizing being tested:
- baseline sizing from the simulation notebook: perfect quality of service (QoS) with zero load shedding
- baseline modified with a halved generator sizing: very good QoS with a bit of load shedding → not penalized
- small PV and small wind generators only: low LCOE (i.e. the production-only LCOE of these sources) but but extremely bad QoS → huge penalty

Documentation for Microgrids.jl
