"""Fit parameters without forecasting, then profile the transmission rate."""
import quantdiffforecast as q
from quantdiffforecast.examples import options_fit_SEIR_flu1918_dist1_3

if __name__ == "__main__":
    opts = options_fit_SEIR_flu1918_dist1_3()
    q.plot_ode_model(opts)                                   # trajectories before fitting
    fit = q.run_fit(opts, 1, 1, 17, seed=1, n_jobs=-1)       # Run_Fit_ODEModel + plotFit_ODEModel
    print(fit.parameters.T, fit.composite, fit.AICcs, sep="\n")

    # Run_ProfileLikelihood_ODEModel(@options_fit_SEIR_flu1918_dist1_3, 1, 17, opts)
    PL = q.run_profile_likelihood(opts, 1, 17, seed=1, paramIndex=[0], alpha=0.05, numstartpoints=3)
    print(PL["referenceStatus"], PL["param"][0]["status"], PL["param"][0]["CI"])
