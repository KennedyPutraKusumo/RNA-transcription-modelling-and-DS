from scipy.stats import f
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np


def total_sum_of_squares(data, std_error=1):
    mean = data.mean()
    tsq = ((data - mean) / std_error) ** 2
    tss = np.sum(tsq)
    return tss


def explained_sum_of_squares(data, predictions, std_error=1):
    mean = data.mean()
    esq = ((predictions - mean) / std_error) ** 2
    ess = np.sum(esq)
    return ess


def RMSE(data, predictions, std_error=1):
    rss = residual_sum_of_squares(data, predictions, std_error)
    rmse = rss / len(data)
    rmse = np.sqrt(rmse)
    return rmse


def residual_sum_of_squares(data, predictions, std_error=1):
    rsq = ((predictions - data) / std_error) ** 2
    rss = np.sum(rsq)
    return rss


def cross_term(data, predictions, std_error=1):
    mean = data.mean()
    first_term = (data - predictions) / std_error
    second_term = (predictions - mean) / std_error
    cross = 2 * np.sum(first_term * second_term)
    return cross


if __name__ == '__main__':
    from S3_calibrate_model import fitting_function, import_data_fit_inputs, import_data_fit_responses

    dataset_selection = "25 November 2022"

    if dataset_selection == "25 November 2022":
        X_numpy = import_data_fit_inputs("2022_11_25_rsm_doe_data.xlsx", "SubsetData")
        rna_data = import_data_fit_responses("2022_11_25_rsm_doe_data.xlsx", "SubsetData")
        mp = [ 1.37815264e-04, 8.34154669e-04,-9.99997778e-01, 6.30866486e-01,
  1.99999978e+00, 7.26133106e-01, 3.26195309e-01, 1.40735543e+00,
  1.00000000e+00]
    elif dataset_selection == "Publication Dataset":
        X_numpy = import_data_fit_inputs("Experimental_Data.xlsx", "ScaledData")
        rna_data = import_data_fit_responses("Experimental_Data.xlsx", "ScaledData")
    else:
        X_numpy = import_data_fit_inputs("2022_11_25_rsm_doe_data.xlsx", "SubsetData")
        rna_data = import_data_fit_responses("2022_11_25_rsm_doe_data.xlsx", "SubsetData")
        mp = [1.09042209e-04, 1.37811096e-03, -7.08819039e-01, 6.20550958e-06,
              1.99999969e+00, 1.11117705e+00, 7.18763807e-01, 1.99999975e+00,
              1.00000000e+00]
    measurement_std_error = 0.1

    rna_predictions = fitting_function(
        X_numpy,
        *mp,
    )

    tss = total_sum_of_squares(rna_data, std_error=measurement_std_error)
    print(f"Total Sum of Squares (TSS): {tss:.2f}")
    ess = explained_sum_of_squares(rna_data, rna_predictions, std_error=measurement_std_error)
    print(f"Explained Sum of Squares (ESS): {ess:.2f}")
    rss = residual_sum_of_squares(rna_data, rna_predictions, std_error=measurement_std_error)
    print(f"Residual Sum of Squares (RSS): {rss:.2f}")
    cross = cross_term(rna_data, rna_predictions, std_error=measurement_std_error)
    print(f"Cross Sum of Squares (CSS): {cross:.2f}")

    rmse = RMSE(rna_data, rna_predictions, std_error=measurement_std_error)
    print(f"Root Mean Squared Error (RMSE): {rmse:.2f}")

    balance_check = tss - ess - rss - cross
    print(f"TSS - ESS - RSS - CSS: {balance_check:.5f}")

    coefficient_of_determination = ess / tss
    print(f"Coefficient of Determination: {coefficient_of_determination:.4f}")
    print("Coefficient of Determination gives indication to goodness/lack of fit, bounded within [0, 1]. The higher the better.")
    coefficient_of_determination2 = 1 - rss / tss
    print(f"Coefficient of Determination: {coefficient_of_determination2:.4f}")
    print("Coefficient of Determination gives indication to goodness/lack of fit, bounded within [0, 1]. The higher the better.")

    print(f"ANOVA Coefficient (ESS / RSS): {ess/rss:.4f}")
    print(f"ANOVA Coefficient gives indication of how much effect the experimental factors have on the RNA concentrations, the higher the coefficiet, the larger the effect.")
    print(f"If the null hypothesis (factors have no effect) is true, statistic will be close to 1.")
    print(f"p-value for the computed value of F-statistic at {len(rna_data)} and {len(rna_predictions)} degrees of freedom: {f.pdf(ess/rss, dfn=len(rna_data), dfd=len(rna_predictions)) * 100:.10f}%")
    print(f"The p-value reported is the probability of obtaining the F-statistic value we obtained purely due to random chance if the null hypothesis is true.")

    fig = plt.figure(figsize=(19, 7))
    axes0 = fig.add_subplot(131)
    axes0.set_title("Data of RNA Yields")
    axes0.scatter(
        np.linspace(0, rna_data.shape[0], rna_data.shape[0]),
        rna_data,
        marker="o",
        c="tab:green",
        alpha=0.5,
        label="Data",
    )
    axes0.plot(
        np.linspace(0, rna_data.shape[0], rna_data.shape[0]),
        rna_predictions,
        marker="1",
        ls="-",
        c="tab:red",
        alpha=0.5,
        label="Predictions",
    )
    axes0.plot(
        [-10, 1000],
        [np.mean(rna_data), np.mean(rna_data)],
        ls="--",
        alpha=0.5,
        c="black",
        label="Mean of Data",
    )
    axes0.set_xlim([-1, rna_data.shape[0] + 1])
    axes0.set_xlabel("Experiment Number")
    axes0.set_ylabel("RNA Yield [g/L]")
    axes0.legend()

    axes1 = fig.add_subplot(132)
    axes1.set_title("Parity Plot")
    axes1.scatter(
        rna_data,
        rna_predictions,
        label="Data Points",
        c="tab:green",
        alpha=0.5,
        marker="o",
    )
    axes1.plot(
        [-1e5, 1e5],
        [-1e5, 1e5],
        ls="--",
        c="black",
        alpha=0.5,
        label="Parity Line",
    )
    axes1.legend()
    axes1.set_xlim([-0.2, 3.2])
    axes1.set_ylim([-0.2, 3.2])
    axes1.set_xlabel("Experimental RNA Yield [g/L]")
    axes1.set_ylabel("Predicted RNA Yield [g/L]")

    axes2 = fig.add_subplot(133)
    axes2.set_title("Residual Plot")
    axes2.scatter(
        rna_predictions,
        (rna_predictions - rna_data),
        label="Data Points",
        c="tab:green",
        alpha=0.5,
        marker="o",
    )
    axes2.plot(
        [-5, 5.0],
        [0, 0],
        c="black",
        ls="--",
        alpha=0.5,
    )
    axes2.set_xlim([-0.2, 3.2])
    axes2.set_xlabel("Predicted RNA Yield [g/L]")
    axes2.set_ylabel("Residuals of RNA Yield [g/L]")

    fig.tight_layout()
    fig.savefig("lack_of_fit.png", dpi=180)

    plt.show()
