# Run of MicrogridsBenchmark.jl on Dell notebook (i7-1165G7)

## System configuration

Results obtained in Ubuntu 22.04, with external power adapter,
using the “Performance” mode in Gnome Settings/Energy.
This setting can have an enourmous influence:
- “Balanced mode”: 180 µs → 260 µs or 180 µs, *it depends!*
- “Power saver” mode: 180 µs → 425 µs

## Results

timing of simulate(mg):  181.459 μs (163 allocations: 759.03 KiB)

detailed timing of simulate(mg):
- operation:  137.699 μs (24 allocations: 753.47 KiB)
- aggregation:  38.350 μs (9 allocations: 144 bytes)
- economics:  2.740 μs (130 allocations: 5.42 KiB)

## Changes

### 2022-10-03: Julia 1.8 update

Rerun performance benchmark with Julia 1.8.2 (from 1.7):
- performance is better (206 µs → 180 µs).
- one less allocation in `economics`, which becomes 3× faster (although it is still, and even more, a small fraction of total time)

### 2022-06-21: Float64 typing

- typing of sum(renewables_production), which fixes the issue with collect(production())
- typing of OperVarsTraj
