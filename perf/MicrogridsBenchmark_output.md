# Run of MicrogridsBenchmark.jl on Dell notebook (i7-1165G7)

2022-06-21 changes: 
* typing of sum(renewables_production), which fixes the issue with collect(production())
* typing of OperVarsTraj

timing of simulate(mg):  206.804 μs (164 allocations: 759.11 KiB)

detailed timing of simulate(mg):
  156.387 μs (24 allocations: 753.53 KiB)
  40.100 μs (9 allocations: 144 bytes)
  7.843 μs (131 allocations: 5.44 KiB)
