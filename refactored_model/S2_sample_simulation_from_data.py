import pandas as pd
import numpy as np



if __name__ == '__main__':
    from S1_model import simulate
    import numpy as np
    data = pd.read_excel(
        "2022_11_25_rsm_doe_data.xlsx",
        "SubsetData",
    )

    C_RNA0 = 1e-3
    exp_idx = 0
    C_NTP0 = 4 * data["NTP [mM]"][exp_idx]
    C_Mg = data["Mg2+ [mM]"][exp_idx]
    C_T7 = data["T7RNAP [mM]"][exp_idx]
    C_DNA = data["DNA template [mM]"][exp_idx]
    tic = np.empty(5)
    tic[:2] = np.array([C_RNA0, C_NTP0])
    tic[2:] = np.array([C_Mg, C_T7, C_DNA])

    spt = np.linspace(0, data["t [hr]"][exp_idx], 51)

    # k_tr = 1e-2
    # k_deg = 1e-2
    # alpha = {
    #     "RNA": -1,
    #     "NTP": 1,
    #     "Mg": 1,
    #     "T7": 1,
    #     "DNA": 1,
    # }
    # beta = {
    #     "RNA": 1,
    #     "Mg": 1,
    # }
    # mp = np.array([
    #     k_tr,
    #     k_deg,
    #     alpha["RNA"],
    #     alpha["NTP"],
    #     alpha["Mg"],
    #     alpha["T7"],
    #     alpha["DNA"],
    #     beta["RNA"],
    #     beta["Mg"],
    # ])
    mp = [ 5.94209500e-05, 2.93645968e-02,-7.77186155e-01, 1.08629793e+00,
  1.06210630e+00, 1.08507241e+00, 8.06923609e-01, 1.07233715e+00,
  1.00000000e+00]

    y = simulate(
        tic,
        spt,
        mp,
    )
    from matplotlib import pyplot as plt
    fig = plt.figure()
    axes = fig.add_subplot(111)
    axes.scatter(
        spt,
        y[:, 0],
        label="RNA Concentration [mM]",
    )
    axes.scatter(
        spt,
        y[:, 1],
        label="NTP Concentration [mM]",
    )
    axes.set_xlabel("Time (hr)")
    axes.set_ylabel("Concentrations [mM]")
    axes.legend()
    fig.tight_layout()
    plt.show()
