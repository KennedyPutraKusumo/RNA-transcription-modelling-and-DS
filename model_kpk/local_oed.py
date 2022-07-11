from pydex.core.designer import Designer
from model import pydex_simulate
import numpy as np


if __name__ == '__main__':
    zero = 1e-5
    mp = np.array([
        4.34,
        5.55e+05,
        1.94e+05,
        1.20e+06,
        zero,
        zero,
        zero,
        zero,
        zero,
    ])

    designer_1 = Designer()
    designer_1.simulate = lambda ti_controls, model_parameters: pydex_simulate(
        ti_controls,
        model_parameters,
        consider_all_params=False,
    )
    designer_1.ti_controls_candidates = designer_1.enumerate_candidates(
        bounds=[
            [0.01, 0.15],     # total Mg conc.    (mol/L)
            [0.01, 0.15],     # initial NTP conc. (mol/L)
            [0.50, 1.50],     # T7RNAP conc.      (mol/L)
        ],
        levels=[
            21,
            21,
            21,
        ],
    )
    # designer_1.model_parameters = mp
    designer_1.model_parameters = mp[[0, 5]]
    designer_1.start_logging()
    designer_1.initialize(verbose=2)

    load_sensitivities = True
    load_atomics = True
    if load_sensitivities:
        # designer_1.load_sensitivity("/local_oed_result/date_2022-3-15/run_2/run_2_atomics_27_cand.pkl")
        # designer_1.load_sensitivity("/local_oed_result/date_2022-3-15/run_1/run_1_atomics_9261_cand.pkl")
        # designer_1.load_sensitivity("/local_oed_result/date_2022-4-8/run_2/run_2_sensitivity_125_cand.pkl")
        # designer_1.load_sensitivity("/local_oed_result/date_2022-4-8/run_1/run_1_sensitivity_125_cand.pkl")
        # designer_1.load_sensitivity("/local_oed_result/date_2022-4-8/run_1/run_1_sensitivity_27_cand.pkl")
        # designer_1.load_sensitivity("/local_oed_result/date_2022-5-3/run_1/run_1_sensitivity_125_cand.pkl")
        designer_1.load_sensitivity("/local_oed_result/date_2022-7-1/run_1/run_1_sensitivity_9261_cand.pkl")
        save_sensitivities = False
    else:
        save_sensitivities = True
    if load_atomics:
        # designer_1.load_atomics("/local_oed_result/date_2022-3-15/run_2/run_2_atomics_27_cand.pkl")
        # designer_1.load_atomics("/local_oed_result/date_2022-3-15/run_1/run_1_atomics_9261_cand.pkl")
        # designer_1.load_atomics("/local_oed_result/date_2022-4-8/run_2/run_2_atomics_125_cand.pkl")
        # designer_1.load_atomics("/local_oed_result/date_2022-4-8/run_1/run_1_atomics_125_cand.pkl")
        # designer_1.load_atomics("/local_oed_result/date_2022-4-8/run_1/run_1_atomics_27_cand.pkl")
        # designer_1.load_atomics("/local_oed_result/date_2022-5-3/run_1/run_1_atomics_125_cand.pkl")
        designer_1.load_atomics("/local_oed_result/date_2022-7-1/run_1/run_1_atomics_9261_cand.pkl")
        save_atomics = False
    else:
        save_atomics = True

    designer_1._norm_sens_by_params = True
    designer_1.design_experiment(
        designer_1.d_opt_criterion,
        package="cvxpy",
        optimizer="MOSEK",
        save_sensitivities=save_sensitivities,
        save_atomics=save_atomics,
        regularize_fim=False,
    )
    fig = designer_1.plot_optimal_controls(non_opt_candidates=False)
    axes = fig.get_axes()[0]
    axes.set_xlim([0.01, 0.15])
    axes.set_ylim([0.01, 0.15])
    axes.set_zlim([0.50, 1.50])
    axes.set_xlabel("Total Mg Concentration (mol/L)")
    axes.set_ylabel("Initial NTP Concentration (mol/L)")
    axes.set_zlabel("T7RNAP Enzyme Concentration (mol/L)")
    # -140 azimuth, 15 elevation
    axes.view_init(15, -140)
    designer_1.print_optimal_candidates()
    designer_1.stop_logging()
    designer_1.show_plots()
