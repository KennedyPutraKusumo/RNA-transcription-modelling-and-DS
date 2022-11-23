from scipy.integrate.odepack import odeint
import numpy as np
import pandas as pd


def dydt(y, t, exp_factors, mp):
    C_RNA = y[0]
    C_NTP = y[1]

    C_Mg = exp_factors[0]
    C_T7 = exp_factors[1]
    C_DNA = exp_factors[2]

    k_tr = mp[0]
    k_deg = mp[1]
    alpha = {
        "RNA": mp[2],
        "NTP": mp[3],
        "Mg": mp[4],
        "T7": mp[5],
        "DNA": mp[6],
    }
    beta = {
        "RNA": mp[7],
        "Mg": mp[8],
    }

    N_all = 1

    V_tr = k_tr * (C_RNA ** alpha["RNA"]) * (C_NTP ** alpha["NTP"]) * (C_Mg ** alpha["Mg"]) * (C_T7 * alpha["T7"]) * (C_DNA ** alpha["DNA"])
    V_deg = k_deg * (C_RNA ** beta["RNA"]) * (C_Mg ** alpha["Mg"])

    dCRNAdt = V_tr - V_deg
    dCNTPdt = -N_all * V_tr
    return np.array([
        dCRNAdt,
        dCNTPdt,
    ])

def simulate(ti_controls, sampling_times, model_parameters):
    y0 = ti_controls[:2]
    exp_factors = ti_controls[2:]
    y = odeint(
        dydt,
        y0=y0,
        t=sampling_times,
        args=(exp_factors, model_parameters),
    )
    return y


if __name__ == '__main__':
    C_RNA0 = 1.00
    C_NTP0 = 7
    y0 = np.array([
        C_RNA0,
        C_NTP0,
    ])

    C_Mg = 4
    C_DNA = 5
    C_T7 = 15
    exp_factors = np.array([
        C_Mg,
        C_DNA,
        C_T7,
    ])
    tic = np.empty(5)
    tic[:2] = y0
    tic[2:] = exp_factors
    spt = np.linspace(0, 10, 101)

    k_tr = 1
    k_deg = 0.1
    alpha = {
        "RNA": -1,
        "NTP": 1,
        "Mg": 1,
        "T7": 1,
        "DNA": 1,
    }
    beta = {
        "RNA": 1,
        "Mg": 1,
    }
    mp = np.array([
        k_tr,
        k_deg,
        alpha["RNA"],
        alpha["NTP"],
        alpha["Mg"],
        alpha["T7"],
        alpha["DNA"],
        beta["RNA"],
        beta["Mg"],
    ])
    mp = [ 1.37815264e-04, 8.34154669e-04,-9.99997778e-01, 6.30866486e-01,
  1.99999978e+00, 7.26133106e-01, 3.26195309e-01, 1.40735543e+00,
  1.00000000e+00]
    y = simulate(tic, spt, mp)

    from matplotlib import pyplot as plt
    fig = plt.figure()
    axes = fig.add_subplot(111)
    axes.scatter(
        spt,
        y[:, 0],
        label=r"$C_{RNA}$"
    )
    axes.scatter(
        spt,
        y[:, 1],
        label="$C_{NTP}$",
    )
    axes.set_xlabel("Time (h)")
    axes.set_ylabel("Concentrations (mM)")
    axes.legend()
    fig.tight_layout()

    plt.show()
