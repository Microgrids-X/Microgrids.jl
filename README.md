# Microgrid.jl

This module simulates the energetic operation of a isolated microgrid, returning economic and operation indicators.

## Installation

To install this package, follow the [instructions for unregistered packages](http://pkgdocs.julialang.org/v1/managing-packages/#Adding-unregistered-packages).

## Microgrid components

- Load
- Diesel generator
- Battery
- Photovoltaic
- Wind turbine

## Operation strategy

- Load-following for battery charging
- Rule-based in the given order:
    1. Renewables (photovoltaic and wind turbine)
    2. Battery
    3. Diesel generator