The `cpic` simulator is able to run 1D particle in cell simulations, with
multiple species.

## Parallel acceleration

The OmpSs-2 framework is used to enable paralellization. Different task are
created, which run in parallel once all dependencies are met.
