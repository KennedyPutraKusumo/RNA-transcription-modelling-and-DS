from S3_calibrate_model import fitting_function, import_data_fit_inputs, import_data_fit_responses
import numpy as np

if __name__ == '__main__':
    dataset_selection = "25 November 2022"
    # dataset_selection = "Publication Dataset"

    if dataset_selection == "25 November 2022":
        X = import_data_fit_inputs("2022_11_25_rsm_doe_data.xlsx", "SubsetData")
        Y = import_data_fit_responses("2022_11_25_rsm_doe_data.xlsx", "SubsetData")
        mp = [ 1.37815264e-04, 8.34154669e-04,-9.99997778e-01, 6.30866486e-01,
  1.99999978e+00, 7.26133106e-01, 3.26195309e-01, 1.40735543e+00,
  1.00000000e+00]
    elif dataset_selection == "Publication Dataset":
        X = import_data_fit_inputs("Experimental_Data.xlsx", "ScaledData")
        Y = import_data_fit_responses("Experimental_Data.xlsx", "ScaledData")
        mp = [ 0.01176829, 0.01459927,-0.54281567, 1.        , 1.0052307 , 1.06075637,
               1.13559156, 1.07186936, 1.        ]
    else:
        X = import_data_fit_inputs("2022_11_25_rsm_doe_data.xlsx", "SubsetData")
        Y = import_data_fit_responses("2022_11_25_rsm_doe_data.xlsx", "SubsetData")
        mp = [1.09042209e-04, 1.37811096e-03, -7.08819039e-01, 6.20550958e-06,
              1.99999969e+00, 1.11117705e+00, 7.18763807e-01, 1.99999975e+00,
              1.00000000e+00]
    Y_pred = fitting_function(X, *mp)

    from matplotlib import pyplot as plt
    fig = plt.figure()
    axes = fig.add_subplot(111)
    axes.scatter(
        Y,
        Y_pred,
        label="Experimental Points",
    )
    y_max = np.max([np.max(Y), np.max(Y_pred)])
    axes.plot(
        [0, y_max],
        [0, y_max],
        label="Parity Line",
        ls="--",
        c="k",
        alpha=0.20,
        lw=2,
    )
    axes.set_xlabel("Experimental RNA Yield [uM]")
    axes.set_ylabel("Predicted RNA Yield [uM]")
    axes.legend()
    fig.tight_layout()
    plt.show()
