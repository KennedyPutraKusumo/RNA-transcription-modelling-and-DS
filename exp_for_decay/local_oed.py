from pydex.core.designer import Designer
from decay_model import pydex_simulate
import numpy as np


if __name__ == '__main__':
    designer_1 = Designer()
    designer_1.simulate = pydex_simulate
    designer_1.ti_controls_candidates = designer_1.enumerate_candidates(
        bounds=[
            [1e-8, 1e-6],       # c_H, i.e., pH of solution
            [1e-3, 3e-3],       #
            [1e0, 3e0],         #
        ],
        levels=[
            11,
            11,
            11,
        ]
    )
    designer_1.model_parameters = np.array([
        1.00000000e+06,
        1.01279115e+00,
        1.00003077e+02,
        8.12962958e-01,
        1.27492412e+00,
    ])
    designer_1.sampling_times_candidates = np.array([
        np.linspace(0, 8, 9)
        for _ in range(designer_1.ti_controls_candidates.shape[0])
    ])
    designer_1.start_logging()
    designer_1.initialize(verbose=2)
    designer_1._norm_sens_by_params = True
    designer_1.design_experiment(
        designer_1.d_opt_criterion,
        optimize_sampling_times=False,
    )
    designer_1.print_optimal_candidates()
    designer_1.stop_logging()
