# ![Microgrids.jl](https://github.com/Microgrids-X/Microgrids-artwork/raw/main/svg/Microgrids-jl.svg)

The `Microgrids.jl` package allows simulating the energetic operation of an isolated microgrid,
returning economic and operation indicators.


## Installation and testing

To install the latest released version of Microgrids.jl: `] add Microgrids`.

To install the lastest development version, follow the [instructions for unregistered packages](http://pkgdocs.julialang.org/v1/managing-packages/#Adding-unregistered-packages).
In particular, you can clone this repository and then, from the Julia command line,
navigate to the folder of your local copy and run:
- either `] add .` to install from the `main` branch of that local repository, but without tracking file changes
- or `] dev .` to install from that local repository *with* tracking file changes (see [dev](https://pkgdocs.julialang.org/v1/managing-packages/#developing) install)

The unit tests which ship with the package can be run with `] test Microgrids`.

## Description of Microgrids.jl

<img alt="Microgrid sizing illustration" src='examples/images/microgrid_h2.png' width="250px">

`Microgrids.jl` can model a microgrid project consisting of:
- One load (described by a time series)
- One dispatchable diesel generator 
- One fuel cell
- One Electrolyzer
- One energy storage (battery)
- Any number of non-dispatchable sources, typically renewable like wind or solar power,
  also modeled from on time series

The energy dispatch at each instant of the simulated operation is a simple
“load following” rule-based control.
The load is power in priority from the dispatchable sources,
then the battery, and only using the dispatchable generator as a last recourse.

`Microgrids.jl` is part of the [Microgrids.X](https://github.com/Microgrids-X/) project
which provides sibling packages in other languages (e.g. in Python)
to better serve the need of different users (e.g. students).

The work on the Julia package specifically focuses on:
- **simulation speed** (about 0.5 ms to evaluate one microgrid project,
  using 1 year of load/solar/wind data at an hourly timestep), way better than pure Python (11 ms for the same task).
  See the [perf](perf) folder for simulation performance benchark.
- **differentiable model**: the exact gradient on the objective function
  (with respect to sizing parameters) can be computed using
  [ForwardDiff](http://www.juliadiff.org/ForwardDiff.jl/).

    - Thanks to our careful treatment of types (and thanks to ForwardDiff and Julia),
      computing the gradient with respect to n=3 parameters with ForwardDiff should run
      faster than with finite differences (which take n+1 = 4× simulation time)

### Documentation

Documentation is provided as [Jupyter](https://jupyter.org/) notebooks within the [examples](examples) folder:
1. Performance simulation in [Microgrid_Wind-Solar.ipynb](examples/Microgrid_Wind-Solar.ipynb):
   1. Describe a Microgrid project using `Microgrids.jl` data structures
   2. Simulate its technical and economic performance
   3. Analyze and display the results
2. Sizing optimization in [Microgrid_sizing_optimization.ipynb](examples/Microgrid_sizing_optimization.ipynb)

### Academic work using Microgrids.jl

First experiments with gradient-based microgrid sizing optimization:

E. de Godoy Antunes, P; Haessig, C. Wang, R. Chouhy Leborgne. “Optimal Microgrid Sizing using Gradient-based Algorithms with Automatic Differentiation”. *ISGT Europe 2022*, Oct 2022, Novi Sad, Serbia. [⟨hal-03370004v2⟩](https://hal.archives-ouvertes.fr/hal-03370004)


## Acknowledgements

The development of Microgrids.jl was first led by Evelise de Godoy Antunes.
She was at that time financed in part by
the Coordenação de Aperfeiçoamento de Pessoal de Nı́vel Superior - Brasil (CAPES) – Finance Code 001,
by Conselho Nacional de Desenvolvimento Cientı́fico e Tecnológico - Brasil (CNPq)
and by the grant “Accélérer le dimensionnement des systèmes énergétiques avec
la différentiation automatique” from [GdR SEEDS (CNRS, France)](https://seeds.cnrs.fr/).
