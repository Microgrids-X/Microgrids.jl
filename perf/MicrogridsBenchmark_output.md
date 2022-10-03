# Run of MicrogridsBenchmark.jl on Dell notebook (i7-1165G7)

## System configuration

Results obtained in Ubuntu 22.04, with external power adapter,
using the “Performance” mode in Gnome Settings/Energy.
This setting can have an enourmous influence:
- “Balanced mode”: 180 µs → 260 µs or 180 µs, *it depends!*
- “Power saver” mode: 180 µs → 425 µs

## Results

timing of simulate(mg):  181.090 μs (165 allocations: 759.12 KiB)

detailed timing of simulate(mg):
- operation:  136.317 μs (24 allocations: 753.47 KiB)
- aggregation:  38.359 μs (9 allocations: 144 bytes)
- economics:  2.940 μs (132 allocations: 5.52 KiB)

## Changes

### 2022-10-03: Component's size type becomes parametrized

Simulation performance is almost unchanged,
with only two extra allocations in economics (130 → 132)

### 2022-10-03: Julia 1.8 update

Rerun performance benchmark with Julia 1.8.2 (from 1.7):
- performance is better (206 µs → 180 µs).
- one less allocation in `economics`, which becomes 3× faster (although it is still, and even more, a small fraction of total time)

### 2022-06-21: Float64 typing

- typing of sum(renewables_production), which fixes the issue with collect(production())
- typing of OperVarsTraj
