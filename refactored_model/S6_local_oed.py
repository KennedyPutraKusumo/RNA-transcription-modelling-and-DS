from S1_model import simulate
from pydex.core.designer import Designer
import numpy as np


if __name__ == '__main__':
    designer = Designer()
    designer.simulate = simulate

    lvls = 5
    designer.ti_controls_candidates = designer.enumerate_candidates(
        bounds=np.array([
            [0.10, 1.00],       # C_RNA0
            [1, 7],             # C_NTP0
            [6, 35],            # C_Mg
            [4, 15],            # C_T7
            [0.03, 0.1],        # C_DNA
        ]),
        levels=np.array([
            1,
            lvls,
            lvls,
            lvls,
            lvls,
        ])
    )
    designer.sampling_times_candidates = np.array([
        np.linspace(0, 2, 11) for _ in designer.ti_controls_candidates
    ])
    designer.model_parameters = np.array([
        1.37815264e-04, 8.34154669e-04, -9.99997778e-01, 6.30866486e-01,
        1.99999978e+00, 7.26133106e-01,  3.26195309e-01, 1.40735543e+00,
        1.00000000e+00
    ])
    designer.start_logging()
    designer.initialize(verbose=2)
    designer.design_experiment(
        designer.d_opt_criterion,
        optimize_sampling_times=True,
        trim_fim=True,
    )
    designer.print_optimal_candidates()
    designer.plot_optimal_predictions(write=True)
    designer.plot_optimal_sensitivities(write=True)
    designer.stop_logging()
    designer.show_plots()
