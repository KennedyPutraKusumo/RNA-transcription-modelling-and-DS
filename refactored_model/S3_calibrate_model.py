from scipy.optimize import curve_fit
from S1_model import simulate
import numpy as np
import pandas as pd


def import_data_fit_responses(file_path, sheet_name):
    data = pd.read_excel(
        file_path,
        sheet_name,
    )
    return data["RNA [uM]"].values

def import_data_fit_inputs(file_path, sheet_name):
    data = pd.read_excel(
        file_path,
        sheet_name,
    )

    X = np.empty((data.shape[0], 6))

    C_NTP0 = 4 * data["NTP [mM]"]
    C_RNA0 = np.ones_like(C_NTP0) * 1e-3
    C_Mg = data["Mg2+ [mM]"]
    C_T7 = data["T7RNAP [mM]"]
    C_DNA = data["DNA template [mg/mL]"]

    X[:, :2] = np.array([C_RNA0, C_NTP0]).T
    X[:, 2:5] = np.array([C_Mg, C_T7, C_DNA]).T
    X[:, -1] = data["t [hr]"]

    return X

def fitting_function(X, k_tr, k_deg, alpha_RNA, alpha_NTP, alpha_Mg, alpha_T7, alpha_DNA, beta_RNA, beta_Mg):
    Y = np.empty(X.shape[0])
    for i, single_X in enumerate(X):

        tic = single_X[:-1]
        spt = np.linspace(0, single_X[-1], 51)
        mp = np.array([
            k_tr,
            k_deg,
            alpha_RNA,
            alpha_NTP,
            alpha_Mg,
            alpha_T7,
            alpha_DNA,
            beta_RNA,
            beta_Mg,
        ])
        y = simulate(
            tic,
            spt,
            mp,
        )
        Y[i] = y[-1, 0]
    return Y


if __name__ == '__main__':
    dataset_selection = "25 November 2022"
    if dataset_selection == "25 November 2022":
        X = import_data_fit_inputs("2022_11_25_rsm_doe_data.xlsx", "SubsetData")
        RNA_data = import_data_fit_responses("2022_11_25_rsm_doe_data.xlsx", "SubsetData")
    elif dataset_selection == "Publication Dataset":
        X = import_data_fit_inputs("Experimental_Data.xlsx", "ScaledData")
        RNA_data = import_data_fit_responses("Experimental_Data.xlsx", "ScaledData")
    else:
        X = import_data_fit_inputs("2022_11_25_rsm_doe_data.xlsx", "SubsetData")
        RNA_data = import_data_fit_responses("2022_11_25_rsm_doe_data.xlsx", "SubsetData")

    k_tr = 1e-2
    k_deg = 1e-2
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
    mp0 = np.array([
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
    mp_bounds = [
        np.array([
            0,
            0,
            -1,
            -1,
            -1,
            -1,
            -1,
            -1,
            -1,
        ]),
        np.array([
            1e-1,
            1e-1,
            -0.1,
            2,
            2,
            2,
            2,
            2,
            2,
        ]),
    ]
    mp_star, cov = calibrate_test = curve_fit(
        fitting_function,
        X,
        RNA_data,
        p0=mp0,
        bounds=mp_bounds,
    )
    print(f"Optimal Parameters".center(100, "="))
    print(np.array2string(calibrate_test[0], separator=","))
    print(f"Covariance".center(100, "="))
    print(np.array2string(calibrate_test[1], separator=","))
