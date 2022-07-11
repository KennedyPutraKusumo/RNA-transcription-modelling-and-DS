from scipy.integrate import odeint
from matplotlib import pyplot as plt
import numpy as np


def dydt(y, t, k, controls):
    # model parameters
    k_ac = k[0]
    n_ac = k[1]
    k_mg = k[2]
    n_mg = k[3]
    n_RNA = k[4]

    # control variables
    c_H = controls[0]
    c_Mg = controls[1]

    # state variable
    c_RNA = y

    V_deg = (k_ac * c_H ** n_ac + k_mg * c_Mg ** n_mg) * c_RNA ** n_RNA

    dRNA_dt = - V_deg

    return dRNA_dt


def model(controls, initial_conditions, sampling_times, model_parameters, plot=False):
    solution = odeint(
        dydt,
        t=sampling_times,
        y0=[initial_conditions, ],
        args=(model_parameters, controls),
    )

    if plot:
        plt.plot(
            t,
            solution[:, 0],
        )
        plt.show()

    c_RNA_predictions = solution[:, 0]

    return c_RNA_predictions


def pydex_simulate(ti_controls, sampling_times, model_parameters):
    controls = [ti_controls[0], ti_controls[1]]
    initial_conditions = ti_controls[2]
    c_RNA = model(
        controls,
        initial_conditions,
        sampling_times,
        model_parameters,
    )
    return c_RNA


if __name__ == '__main__':

    # controls
    c_H = 1e-7      # pH of 7: pH = -log10([H+]) -> [H+] = 10^(-pH)
    c_Mg = 1e-3     # total Mg conc. taken from original DS publication
    controls = np.array([
        c_H,
        c_Mg,
    ])

    # initial conditions
    c_RNA_0 = 1e0

    # model parameters
    mp = np.array([
        1e6,        # k_ac
        1,          # n_ac
        1e2,        # k_mg
        1,          # n_mg
        1,          # n_RNA
    ])

    # sampling times
    t = np.linspace(0, 30, 11)

    c_RNA_prediction = model(
        controls,
        c_RNA_0,
        t,
        mp,
        plot=True,
    )

