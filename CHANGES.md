# CHANGELOG for Microgrids.jl

## v0.11.0

Jul 30, 2024

Cumulated changes of Spring/Summer 2024:

- Implement **economically consistent salvage value** definition, see our report Pierre Haessig. Economic consistency of salvage value definitions. 2024. [⟨hal-04097092v2⟩](https://hal.science/hal-04097092v2)
  - Microgrid's `Project` description structure accepts a new `salvage_type` argument, with the default value being `LinearSalvage` which is the classical and previous behavior, but which can be changed to `ConsistentSalvage`
- Add new `smoothing` parameter of the toplevel microgrid `simulate` function. It replaces and generalizes the relaxation parameter `ε` introduced in v0.10.2. 
  - Default value is `NoSmoothing` but can be set to any `Smoothing` instance
  - the former`ε` parameter corresponds to the first field of the `Smoothing` structure, which is named `transition`. 
  - the 2nd parameter if the `Smoothing` structure is `gain`, which is 1.0 by default, but can be increased to heuristically compensate for the underapproximation when using a strictly positive `transition` value.
- The structures for describing Microgrid projects (`Microgrid`, `Project` and all components such as `Battery`) are now:
  - **mutable** (fields can be changed after creation)
  - accept **keyword** arguments
  - have sensible **defaults** for secondary parameters, matching Microgrids.py interface 


## v0.10.2

Oct 13, 2023

- add relaxation parameter `ε` to toplevel microgrid `simulate` function.
  - this allow reproducing the optimization results from ou article: de Godoy Antunes, P. Haessig, C. Wang, and R. Chouhy Leborgne, “Optimal  Microgrid Sizing using Gradient-based Algorithms with Automatic  Differentiation,” ISGT Europe 2022, Novi Sad, Serbia, 2022. https://dx.doi.org/10.1109/ISGT-Europe54678.2022.9960498, Archived at https://hal.archives-ouvertes.fr/hal-03370004
  - but a reimplementation of these results in Microgrids.jl example folder is yet to be done.

## v0.10.1

Sep 27, 2023

Same as [v0.10.0](https://github.com/Microgrids-X/Microgrids.jl/releases/tag/v0.10.0) + minor fixes to performance benchmark and README

## v0.10.0

Sep 27, 2023

Microgrids.jl updated to:

- match Microgrids.py interface (same field names in data structure)
- add wind power model
- improved notebook examples, with both simulation and optimization

## v0.9.0: first registered version of Microgrids.jl

Nov 7, 2022

Now Microgrids.jl should be installable with `] add Microgrids`!
