## A postprocessing script to compare Batree simulations with PyBAMM

Install PyBAMM with `pip install pybamm` and Matplotlib with `pip install matplotlib`.

Choose what to plot by modifying script variables:
```
PLOT_SPM   = True
PLOT_SPMe  = True
PLOT_P2D   = True

PLOT_MFEM  = True
PLOT_PYBAMM = True
```

Ensure variable defining the mfem executable is correct:
```
mfem_executable = "./../batree"
```

Note: Batree simulations seem to take quite a long time (approx 5 mins).  Not currently tested in parallel.
