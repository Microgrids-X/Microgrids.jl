# Run of MicrogridsBenchmark.jl on Dell notebook (i7-1165G7)

## System configuration

Julia version 1.10.4

Results obtained in Ubuntu 22.04, with external power adapter,
using the “Performance” mode in Gnome Settings/Energy.
This setting can have an major influence:
- “Balanced” mode: often yields the best performance, but sometimes quite variable
- “Power saver” mode: more than twice slower (144 µs → 344 µs)

## Results

timing of simulate(mg):  145 to 148 μs (122 allocations: 689.00 KiB)

detailed timing of simulate(mg):
- operation:  120    μs (22 allocations: 684.98 KiB)
- aggregation: 21  μs (19 allocations: 304 bytes)
- economics:    1.9 μs (81 allocations: 3.72 KiB)

timing of gradient(sim_npc, x):  466 to 480 μs (128 allocations: 2.68 MiB)
which represents 3.2× simulation time (gradient of dim 3).

## Changes

### 2024-07-30: small performance decrease due to warm environment?

Test run in room at 24°C

### 2024-06-11: small performance decrease due to ConsistentSalvage

3 more allocations in economics, but economics is still < 2µs

### 2023-07-25: small performance increase due to refactor

About 8 µs saved in the refactor since 2023-07-18:
- operation: 6 µs save due to a one line change?
- economics: twice as fast

### 2023-07-18: Updgrade Julia 1.8 → 1.9

Small performance *decrease* with the Julia upgrade 1.8 → 1.9
on `simulate(mg)`: 144 → 153 µs (+6%),
with most increase coming from `operation`.
Still, in relative terms, `economics` also slowed down: 2.9 → 3.3 µs

Stronger negative impact on gradient computation: 390 µs → 550 µs (+41%),
which represents now 3.6× the simulation time, almost as slow as finite difference.
(remark: ForwardDiff was also updated in the process)

### 2023-07-14: Refactor operation and aggregation

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
