import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import pybamm
import subprocess

# Close all figures
plt.close("all")

# Plotting directory, make if it doesn't exist
plots_dir = "./plots/"
os.makedirs(plots_dir, exist_ok=True)

mfem_executable = "./../batree"

# Plot flags
PLOT_SPM   = True
PLOT_SPMe  = True
PLOT_DFN   = True

PLOT_MFEM  = True
PLOT_PYBAMM = True

# Set parameter values for PyBaMM simulation (Chen2020 is the LG M50 battery)
parameter_values = pybamm.ParameterValues("Chen2020")

# Colours
red   = (0.7, 0.2, 0.2)
blue  = (0.2, 0.2, 0.7)
black = (0.2, 0.2, 0.2)

# For reducing number of markers on plot.
mfem_indx = 30

# Standard plot settings for all plots.
def easyplot(x, y, colour, linestyle, fig, label, marker_indx):
    plt.figure(fig)
    plt.plot(
        x, y,
        linestyle,
        color=colour,
        linewidth=1.3,
        label=label,
        markerfacecolor="w",
        markevery=marker_indx
    )

# Function to run Batree and extract results from output.
def run_batree(mfem_executable,sim_type):

    print(''.join(["Running ", sim_type, " simulation in batree..."]))
    cmd = [mfem_executable, "-m", sim_type]
    result = subprocess.run(cmd, capture_output=True, text=True)

    times = []
    voltages = []

    for line in result.stdout.splitlines():
        parts = line.split()

        if len(parts) == 3 and parts[0].isdigit():
            step, time, voltage = parts
            times.append(float(time))
            voltages.append(float(voltage))

    return times, voltages


# Function to run PyBAMM.
def run_pybamm(model, parameter_values):

    print(''.join(["Running ", str(model.__class__.__name__), " simulation in PyBaMM..."]))
    sim = pybamm.Simulation(model, parameter_values=parameter_values)
    soln = sim.solve([0, 3600])

    time = soln["Time [s]"].entries
    voltage = soln["Voltage [V]"].entries

    return time, voltage


# Running simulations and plotting results.
if PLOT_SPM and PLOT_MFEM:

    time, voltage = run_batree(mfem_executable, "SPM")
    easyplot(time, voltage, blue, "-o", 1, "SPM (MFEM)", mfem_indx)

if PLOT_SPMe and PLOT_MFEM:

    time, voltage = run_batree(mfem_executable, "SPMe")
    easyplot(time, voltage, red, "-o", 1, "SPMe (MFEM)", mfem_indx)

if PLOT_DFN and PLOT_MFEM:

    time, voltage = run_batree(mfem_executable, "P2D")
    easyplot(time, voltage, black, "-o", 1, "DFN (MFEM)", mfem_indx)

if PLOT_SPM and PLOT_PYBAMM:

    time,voltage = run_pybamm(pybamm.lithium_ion.SPM(), parameter_values)
    easyplot(time, voltage, blue, "--s", 1, "SPM (PyBaMM)", 2)

if PLOT_SPMe and PLOT_PYBAMM:

    time,voltage = run_pybamm(pybamm.lithium_ion.SPMe(), parameter_values)
    easyplot(time, voltage, red, "--s", 1, "SPMe (PyBaMM)", 5)

if PLOT_DFN and PLOT_PYBAMM:

    time,voltage = run_pybamm(pybamm.lithium_ion.DFN(), parameter_values)
    easyplot(time, voltage, black, "--s", 1, "DFN (PyBaMM)", 5)


# Tidying up and saving figure.
plt.figure(1)
plt.xlim(0, 3600)
plt.xlabel("Time (s)")
plt.ylabel("Voltage (V)")
plt.legend(loc="lower left")
plt.savefig(os.path.join(plots_dir, "voltage_compare.png"), dpi=300)
