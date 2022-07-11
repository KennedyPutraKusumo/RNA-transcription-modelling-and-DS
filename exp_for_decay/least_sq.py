from scipy.optimize import minimize
from decay_model import model
from matplotlib import pyplot as plt
import numpy as np


def lsq(mp, data):
    controls = data["controls"]
    initial_conditions = data["initial_conditions"]
    sampling_times = data["sampling_times"]
    c_RNA_measurements = data["c_RNA_measurements"]

    residuals = []
    for experiment in zip(controls, initial_conditions, sampling_times,
                          c_RNA_measurements):
        predictions = model(
            experiment[0],
            experiment[1],
            experiment[2],
            mp,
        )
        measurement = np.array(experiment[3])
        residual = predictions - measurement
        residuals.append(residual)
    residuals = np.array(residuals)
    sq_residuals = np.square(residuals)
    ssq = np.sum(sq_residuals)

    return ssq


def calibrate(mp0, data, bounds=None):
    opt_result = minimize(
        lsq,
        mp0,
        args=data,
        # method="l-bfgs-b",
        method="SLSQP",
        bounds=bounds,
        options={
            "disp": True,
        }
    )
    calibrated_mp = opt_result.x
    return calibrated_mp


if __name__ == '__main__':
    half_life = 1.5     # hours
    artificial_data = {
        "controls": [
            [1e-7, 1e-3],
            [1e-7, 1e-3],
        ],
        "initial_conditions": [
            1e0,
            2e0,
        ],
        "sampling_times": [
            [0, half_life, half_life * 2, half_life * 3, half_life * 4],
            [0, half_life, half_life * 2, half_life * 3, half_life * 4],
        ],
        "c_RNA_measurements": [
            [1.0, 0.50, 0.25, 0.125, 0.0625],
            [2.0, 1.00, 0.50, 0.250, 0.1250],
        ],
    }

    mp = np.array([
        1e6,        # k_ac
        1,          # n_ac
        1e2,        # k_mg
        1,          # n_mg
        1,          # n_RNA
    ])

    res = lsq(
        mp,
        artificial_data,
    )

    estimated_mp = calibrate(
        mp,
        artificial_data,
        bounds=[
            (1e5, 1e8),
            (1, 2),
            (1e1, 1e4),
            (0.1, 4),
            (0.1, 4),
        ],
    )

    fig = plt.figure()

    axes1 = fig.add_subplot(111)

    axes1.scatter(
        artificial_data["sampling_times"][0],
        artificial_data["c_RNA_measurements"][0],
        c="tab:red",
        label="Experiment 1"
    )
    axes1.scatter(
        artificial_data["sampling_times"][1],
        artificial_data["c_RNA_measurements"][1],
        c="tab:blue",
        label="Experiment 2"
    )

    initial_predictions = model(
        artificial_data["controls"][0],
        artificial_data["initial_conditions"][0],
        artificial_data["sampling_times"][0],
        mp,
    )
    axes1.plot(
        artificial_data["sampling_times"][0],
        initial_predictions,
        c="tab:red",
        alpha=0.25,
        ls="--",
        label="Before Calibration",
    )

    calibrated_predictions = model(
        artificial_data["controls"][0],
        artificial_data["initial_conditions"][0],
        artificial_data["sampling_times"][0],
        estimated_mp,
    )
    axes1.plot(
        artificial_data["sampling_times"][0],
        calibrated_predictions,
        c="tab:red",
        label="After Calibration",
    )

    initial_predictions = model(
        artificial_data["controls"][1],
        artificial_data["initial_conditions"][1],
        artificial_data["sampling_times"][1],
        mp,
    )
    axes1.plot(
        artificial_data["sampling_times"][1],
        initial_predictions,
        c="tab:blue",
        alpha=0.25,
        ls="--",
        label="Before Calibration",
    )

    calibrated_predictions = model(
        artificial_data["controls"][1],
        artificial_data["initial_conditions"][1],
        artificial_data["sampling_times"][1],
        estimated_mp,
    )
    axes1.plot(
        artificial_data["sampling_times"][1],
        calibrated_predictions,
        c="tab:blue",
        label="After Calibration",
    )
    axes1.legend()
    axes1.set_xlabel("Time (hours)")
    axes1.set_ylabel("[RNA] (mol/L)")

    print(f"Estimated parameters are: {estimated_mp}")

    plt.show()
