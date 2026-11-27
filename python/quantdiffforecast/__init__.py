"""QuantDiffForecast for Python.

A port of the MATLAB toolbox QuantDiffForecast (G. Chowell, A. Bleichrodt,
R. Luo; Statistics in Medicine 2024, doi:10.1002/sim.10036): parameter
estimation and forecasting with quantified uncertainty for ODE models.

Main entry points (MATLAB name -> Python name)

    Run_Fit_ODEModel                  -> run_fit
    Run_Forecasting_ODEModel          -> run_forecasting
    plotFit_ODEModel / plotForecast_ODEModel -> plot_fit / plot_forecast
    plotODEModel                      -> plot_ode_model
    Run_ProfileLikelihood_ODEModel    -> run_profile_likelihood
    quantdiffCompareForecasts         -> compare_forecasts
    run_ensemble                      -> run_ensemble
    baselineForecast                  -> baseline_forecast
    quantdiffForecastScores           -> forecast_scores
    Run_simulate_ODEModel             -> run_simulate
    Run_PracticalIdentifiability_ODEModel -> run_practical_identifiability
    summarizeIdentifiability          -> summarize_identifiability
"""
from .baseline import baseline_forecast
from .compare import compare_forecasts, ensemble_weights, linear_pool, run_ensemble
from .config import Model, Options, Params, QuantDiffForecastError, Vars
from .core import (ObjectiveContext, aicc_count_adjustment, effective_n, get_aicc, get_nparams,
                   initial_conditions, objective_value, observation_curve, ode_tolerances, solve_ode)
from .identifiability import bound_hit_fractions, pii, quantile_mcse
from .noise import (add_error_structure, estimate_laplace_scale, estimate_pearson_dispersion,
                    negative_noise_share)
from .optimizer import FitResult, fit_model
from .plotting import (plot_fit, plot_forecast, plot_ode_model, plot_parameters, plot_profiles,
                       plot_states)
from .profile import run_profile_likelihood
from .runners import RunResult, load_forecast, run_fit, run_forecasting
from .scoring import (compute_forecast_performance, compute_quantiles, compute_wis,
                      forecast_scores, mquantile, scores_from_quantiles)
from .simulate import run_practical_identifiability, run_simulate, summarize_identifiability
from . import examples, models

__version__ = "0.1.1"
