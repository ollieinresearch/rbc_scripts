# Python API

Every module of the `rbc` package, generated from the docstrings. Each module is also a
command where it says so; its options are in {doc}`cli`. Modules needing MPI, Dedalus or
VTK are documented without importing those packages, so these pages build anywhere.

```{eval-rst}
.. autosummary::
   :toctree: api/generated
   :recursive:

   rbc
```

## The simulation script

`sim/rayleigh_benard.py` is not part of the package, but its helpers are documented here.

```{eval-rst}
.. automodule:: rayleigh_benard
   :members: get_parser, log_job, main
```
