from odes_and_curve_fitting_functions import datafitting_transcription_experimental
from matplotlib import pyplot as plt
import numpy as np


def pydex_simulate(ti_controls, model_parameters, consider_all_params=False):
    X = np.array([[
        6,                  # batch time fixed to 6 hours
        ti_controls[0],     # total Mg conc.    (mol/L)
        ti_controls[1],     # initial NTP conc. (mol/L)
        2e-4,               # spermidine concentration fixed to 2E-4
        ti_controls[2],     # T7RNAP conc.      (mol/L)
    ]])
    # =============================
    # consider all model parameters
    # =============================
    if consider_all_params:
        rna_yields = datafitting_transcription_experimental(
            X=X,
            k_app=model_parameters[0],
            K1=model_parameters[1],
            K2=model_parameters[2],
            k_ac=model_parameters[3],
            k_ba=model_parameters[4],
            k_Mg=model_parameters[5],
            K3=model_parameters[6],
            K4=model_parameters[7],
            K5=model_parameters[8],
        )
    # =============================
    # consider only identifiable: k_app, k_mg
    # =============================
    else:
        rna_yields = datafitting_transcription_experimental(
            X=X,
            k_app=model_parameters[0],
            K1=5.55e+05,
            K2=1.94e+05,
            k_ac=1.20e+06,
            k_ba=0,
            k_Mg=model_parameters[1],
            K3=0,
            K4=0,
            K5=0,
        )
    return rna_yields[0, 0, :]


if __name__ == '__main__':
    mp = np.array([
        4.34,
        5.55e+05,
        1.94e+05,
        0,
        0,
        0,
        0,
        0,
        0,
    ])

    # test simulate function: single simulation
    if True:
        tic = np.array([
            0.01,
            0.01,
            0.50,
        ])

        y = pydex_simulate(
            tic,
            mp,
        )
        print(y)
    # test simulate function: simulate yields over grids of (i) x1, x2 (ii) x1, x3
    if True:
        grid_reso = 11j
        x1, x2 = np.mgrid[0.01:0.15:grid_reso, 0.01:0.15:grid_reso]
        x1 = x1.flatten()
        x2 = x2.flatten()
        x3 = np.ones_like(x1) * 0.50
        tics = np.array([x1, x2, x3]).T
        yields = []
        for tic in tics:
            yields.append(pydex_simulate(
                tic,
                mp,
                consider_all_params=True,
            ))
        yields = np.array(yields)
        fig = plt.figure(figsize=(12, 5))
        axes1 = fig.add_subplot(121, projection="3d")
        axes1.scatter(
            tics[:, 0],
            tics[:, 1],
            yields[:, 0],
        )
        axes1.set_xlabel("Total Mg conc.(mol/L)")
        axes1.set_ylabel("Initial NTP conc. (mol/L)")
        axes1.set_zlabel("Effective RNA yield (g/L)")

        x1, x3 = np.mgrid[0.01:0.15:grid_reso, 0.50:1.50:grid_reso]
        x1 = x1.flatten()
        x2 = np.ones_like(x1) * 0.01
        x3 = x3.flatten()
        tics = np.array([x1, x2, x3]).T
        yields = []
        for tic in tics:
            yields.append(pydex_simulate(
                tic,
                mp,
                consider_all_params=True,
            ))
        yields = np.array(yields)
        axes2 = fig.add_subplot(122, projection="3d")
        axes2.scatter(
            tics[:, 0],
            tics[:, 2],
            yields[:, 0],
        )
        axes2.set_xlabel("Total Mg conc.(mol/L)")
        axes2.set_ylabel("T7RNAP conc. (mol/L)")
        axes2.set_zlabel("Effective RNA yield (g/L)")
        fig.tight_layout()
        plt.show()
