from odes_and_curve_fitting_functions import datafitting_transcription_experimental
import numpy as np


def log_lkhd(theta, x, y, y_err):
    model = datafitting_transcription_experimental(
        x,
        *theta,
        k_ac=0,
        k_ba=0,
        k_Mg=0,
        K3=0,
        K4=0,
        K5=0,
    )
    return np.sum(((y - model[0, :, 0])/y_err) ** 2)

def log_prior(theta):
    if 0 <= theta[0] <= 5 and 1e5 <= theta[1] <= 1e6 and 1e5 <= theta[2] <= 1e6:
        return 0
    else:
        return -np.inf

def log_probability(theta, x, y, y_err):
    prior = log_prior(theta)
    if not np.isfinite(prior):
        return -np.inf
    else:
        return prior + log_lkhd(theta, x, y, y_err)


if __name__ == '__main__':
    from matplotlib import pyplot as plt
    import pandas as pd
    import emcee as mc

    lsq_mp = [4.34, 5.55e+05, 1.94e+05]
    pos = lsq_mp + 1e-4 * np.random.randn(32, 3)
    nwalkers, ndim = pos.shape

    exp_data = pd.read_excel(
        "Experimental_Data.xlsx",
        sheet_name="ModelData",
        header=0,
    )
    print(list(exp_data))
    print(exp_data.head())

    rna_data = exp_data["RNA [g/L]"]
    measurement_std_error = 0.1

    X = exp_data.loc[:, "t [hr]":"T7RNAP ratio"]
    X_numpy = X.to_numpy()
    sampler = mc.EnsembleSampler(
        nwalkers,
        ndim,
        log_probability,
        args=(X_numpy, rna_data.to_numpy(), measurement_std_error)
    )
    sampler.run_mcmc(pos, 100, progress=True)

    fig, axes = plt.subplots(3, figsize=(10, 7), sharex=True)
    samples = sampler.get_chain()
    labels = ["m", "b", "log(f)"]
    for i in range(ndim):
        ax = axes[i]
        ax.plot(samples[:, :, i], "k", alpha=0.3)
        ax.set_xlim(0, len(samples))
        ax.set_ylabel(labels[i])
        ax.yaxis.set_label_coords(-0.1, 0.5)

    axes[-1].set_xlabel("step number")

    tau = sampler.get_autocorr_time()
    print(tau)

    flat_samples = sampler.get_chain(discard=100, thin=15, flat=True)
    print(flat_samples.shape)

    import corner

    fig = corner.corner(
        flat_samples, labels=labels, truths=lsq_mp,
    )

    inds = np.random.randint(len(flat_samples), size=100)
    for ind in inds:
        sample = flat_samples[ind]
        plt.plot(x0, np.dot(np.vander(x0, 2), sample[:2]), "C1", alpha=0.1)
    plt.errorbar(x, y, yerr=yerr, fmt=".k", capsize=0)
    plt.plot(x0, m_true * x0 + b_true, "k", label="truth")
    plt.legend(fontsize=14)
    plt.xlim(0, 10)
    plt.xlabel("x")
    plt.ylabel("y")

    plt.show()
