# Run of MicrogridsBenchmark.jl

## System configuration

- Computer: Dell notebook i7-1165G7
- Julia version 1.9.2

Results obtained in Ubuntu 22.04, with external power adapter,
using the “Performance” mode in Gnome Settings/Energy.
This setting can have an major influence:
- “Balanced” mode: often yields the best performance, but sometimes quite variable
- “Power saver” mode: more than twice slower (144 µs → 344 µs)

## Results

timing of simulate(mg):  90 μs (61 allocations: 71.86 KiB)

detailed timing of simulate(mg):
- operation: 89   μs (5 allocations: 68.72 KiB)
- economics: 1.37 μs (56 allocations: 3.14 KiB)

timing of gradient(sim_npc, x): 154 μs (67 allocations: 279.78 KiB)
which represents 1.7× simulation time (gradient of dim 3).

## Changes

### 2023-07-24: Combined operation-aggregation

Operation now includes aggregation (like in Microgrids.py) which saves allocating
several operational time series which are then to be aggregated.

Result:
- good simulation speed up (155 → 90 µs)
- much smaller memory allocation (691 → 72 KiB)
- even better speed up for differentiated simulation (550 → 154 μs)

Remark: there is still one remaining source of vector allocations
for computing the production of each NonDispatchable source.
See vectorless-operation branch for the work in progress to remove these last vectors.

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
