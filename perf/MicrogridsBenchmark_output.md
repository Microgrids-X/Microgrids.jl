# Run of MicrogridsBenchmark.jl on Dell notebook (i7-1165G7)

**Benchmark results outdated** -> to be run again on Dell i7

## System configuration

Julia version 1.8.3

Results obtained in Ubuntu 22.04, with external power adapter,
using the “Performance” mode in Gnome Settings/Energy.
This setting can have an major influence:
- “Balanced” mode: often yields the best performance, but sometimes quite variable
- “Power saver” mode: more than twice slower (144 µs → 344 µs)

## Results

timing of simulate(mg):  144 to 146 μs (173 allocations: 691.25 KiB)

detailed timing of simulate(mg):
- operation:  117 μs (22 allocations: 684.98 KiB)
- aggregation: 20.3 μs (19 allocations: 304 bytes)
- economics:  2.93 μs (130 allocations: 5.77 KiB)

timing of gradient(sim_npc, x):  388 μs (187 allocations: 2.69 MiB),
which represents 2.7× simulation time (gradient of dim 3).

## Changes

### 2022-07-14: Refactor operation and aggregation

Good news: about 40 µs saved on simulate(mg): 181 → 144 µs
- about 20 µs saved on operation: 136 → 117 µs
  (perhaps thanks to the 68 KiB allocation saved, i.e. one hourly time series?)
- about 20 µs saved on aggregation: 38 → 20 µs

Perhaps even faster simulation could be achievied by merging aggregation
into operation (i.e. make time series recording optional), like in Microgrids.py.

Not so good: gradient evaluation timing increased from 353 → 388 µs,
albeit with quite some variability (despite ~60 KiB allocation saved),
so the timing ration gradient/simulation increases to 2.7.

### 2022-10-25: Add gradient timing

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
