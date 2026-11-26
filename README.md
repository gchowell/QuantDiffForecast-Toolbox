# QuantDiffForecast

QuantDiffForecast fits ordinary differential equation (ODE) models to one or more observed time series in MATLAB. It estimates parameters, generates short-term forecasts with prediction intervals, and evaluates forecasts against held-out observations. The quick start uses daily influenza counts from San Francisco in 1918. A second example fits reported COVID-19 cases in Spain. You can use your own model by specifying its equations and how each observed series relates to a model state.

[Quick start](#quick-start) · [Use your own data](#using-your-own-data) · [Forecast evaluation](#read-a-forecast-score) · [Compare forecasts](#compare-models-across-forecast-origins) · [Initial conditions](#initial-conditions-for-each-window) · [Profile likelihood](#profile-likelihood) · [Parameter recovery](#practical-identifiability-simulations) · [Output files](#outputs) · [Figures](#example-visualizations) · [Published tutorial](https://doi.org/10.1002/sim.10036) · [Video tutorial](https://www.youtube.com/watch?v=eyyX63H12sY&t=41s)

## Requirements and installation

The MATLAB source requires **Optimization Toolbox** (`fmincon`), **Global Optimization Toolbox** (`MultiStart`), and **Statistics and Machine Learning Toolbox** (probability distributions and random sampling). **Parallel Computing Toolbox** is optional for parallel optimization and simulation studies.

Download the ZIP from the repository’s **Code** menu, or clone it in a terminal.

```bash
git clone https://github.com/gchowell/QuantDiffForecast-Toolbox.git
```

In MATLAB, set **Current Folder** to the repository root and run the following. Examples below assume the code folder is your current folder. The runners create `output/` there.

```matlab
codeDir = fullfile(pwd, 'forecasting_odemodels code');
assert(isfolder(codeDir), 'Set Current Folder to the repository root first.');
addpath(codeDir);
cd(codeDir);
```

For a short installation check, run `qdfSmokeTest(2)`. It fits and plots a forecast using 20 bootstrap datasets and prints a verdict. The full example below uses 300. Separate Windows and Mac packages are available. The Windows executable requires MATLAB Runtime R2023b. [Smoke test](forecasting_odemodels%20code/qdfSmokeTest.m) · [Standalone packages](stand%20alone%20executable/) · [Windows deployment instructions](stand%20alone%20executable/readme.txt)

## Quick start

### Forecast ten days from the first 17 observations

The bundled susceptible–exposed–infectious–recovered (SEIR) model fits daily counts to increments of state 5, which accumulates transitions from exposed to infectious. It estimates the transmission rate β and negative-binomial dispersion α. The progression rate κ, recovery rate γ and population size N are fixed. The first observed count sets the initial value of state 5. [Data](forecasting_odemodels%20code/input/curve-flu1918SF.txt) · [Model equations](forecasting_odemodels%20code/SEIR1.m) · [Forecast options](forecasting_odemodels%20code/options_forecast_SEIR_flu1918_dist1_3.m)

```matlab
forecastOptions = @options_forecast_SEIR_flu1918_dist1_3;
rng(1, 'twister');
Run_Forecasting_ODEModel(forecastOptions, 1, 1, 17, 10);
```

The two `1` arguments select one window starting at row 1. The runner fits the observed data, refits 300 bootstrap datasets and generates 20 predictive draws per refit. The resulting 6,000 draws include parameter and observation uncertainty.

| Period | MATLAB rows | Time labels in the file | Observed counts at the endpoints |
| --- | --- | --- | --- |
| Calibration | 1–17 | 0–16 | 4 and 58 |
| Held-out forecast targets | 18–27 | 17–26 | 124 and 1,104 |

Read the predictions alongside the score for day 10 alone and the average score over all ten held-out days.

```matlab
suffix = ['model_name-SEIR model-fit_index-5-fixI0-1-method-3-dist-3-' ...
          'tstart-1-tend-1-calibrationperiod-17-horizon-10-1918 Flu-cases.csv'];
trajectory = readtable(fullfile('output', ['Forecast-' suffix]));
scores = readtable(fullfile('output', ['performance-forecasting-' suffix]));
byLead = readtable(fullfile('output', ['scores-byhorizon-' suffix]));
disp(trajectory(18:27, {'time','data','median','LB','UB','model_median'}))
disp(byLead(end, {'horizon','observed','median','abs_error','cover95','WIS'}))  % day 10
disp(scores(end, :))   % days 1–10
```

`median`, `LB` and `UB` summarize predicted observations (median and central 95% prediction interval). `model_median` summarizes the refitted model trajectories before observation noise is added. To view the figures, run `plotForecast_ODEModel(forecastOptions, 1, 1, 17, 10)` **after reading the CSV**. The plotter rewrites that file with the first five columns. The six-column table above remains available as `trajectory` in MATLAB. [Forecast runner](forecasting_odemodels%20code/Run_Forecasting_ODEModel.m) · [Forecast plotter](forecasting_odemodels%20code/plotForecast_ODEModel.m)

### Fit parameters without forecasting

Use the fitting runner when the analysis ends at the last calibration observation. It accepts a fitting-options function. The forecast example above performs its own fit.

```matlab
fitOptions = @options_fit_SEIR_flu1918_dist1_3;
rng(1, 'twister');
Run_Fit_ODEModel(fitOptions, 1, 1, 17);
plotFit_ODEModel(fitOptions, 1, 1, 17);
```

To inspect trajectories before fitting, use `plotODEModel(fitOptions)`. Both named quick-start options use method 3, 300 bootstrap datasets and `numstartpoints = 10`. Their starting values and lower bounds for β differ. Calling a runner **without arguments** selects a different options file. The default fit uses Poisson likelihood and 300 bootstrap datasets. The default forecast uses method 3 and 100 bootstrap datasets. [Named fit options](forecasting_odemodels%20code/options_fit_SEIR_flu1918_dist1_3.m) · [Default fit](forecasting_odemodels%20code/options_fit.m) · [Default forecast](forecasting_odemodels%20code/options_forecast.m)

The Spain example fits reported cases using an SEIR model with an unreported infectious compartment. To fit the first 60 observations with a negative-binomial likelihood, run `Run_Fit_ODEModel(@options_fit_SEIR_unreported_covid_spain_dist1_3, 1, 1, 60)`. A Poisson fit is also provided. [Spain data](forecasting_odemodels%20code/input/curve-cases-covid-spain.txt) · [Negative-binomial options](forecasting_odemodels%20code/options_fit_SEIR_unreported_covid_spain_dist1_3.m) · [Poisson options](forecasting_odemodels%20code/options_fit_SEIR_unreported_covid_spain_dist1_1.m)

## Using your own data

1. Place a numeric, header-free `.txt` file in the [input folder](forecasting_odemodels%20code/input/). The first column is a consecutive, unit-spaced time index. Each remaining column holds one observed series. Use the observation interval as the time unit for the ODE rates. Resolve missing or irregular observations before fitting.
2. Copy the fit or forecast options file from the quick start. Rename both the file and its first function declaration, keep the output list, and set `cadfilename1` to the input filename. The runners append `.txt` if omitted. Give the analysis a distinct `model.name` so its exports are easy to identify.
3. Map the observed columns to model states, in the same order as the file. Set the model, bounds and observation model in [Configuration](#configuration), then pass the copied options function to the corresponding runner.

### Match observations to model states

For example, suppose the input has columns `time`, `infectious prevalence` and `new cases`. If state 3 is current prevalence and state 5 accumulates new cases, use the following mapping.

```matlab
vars.fit_index = [3, 5];
vars.fit_diff  = [0, 1];
```

| Observed series | Flag | What the toolbox compares with the data |
| --- | ---: | --- |
| State level, such as prevalence or cumulative cases | `0` | The state value at each observation time. |
| Count per interval, such as daily new cases | `1` | `abs([C(1); diff(C)])` for the selected state `C`. The first value uses the initial state. Later values use absolute increments. |

For daily incidence, supply daily counts and choose a nondecreasing cumulative model state. For cumulative observations, use a level mapping. The runner reads each observed column as supplied. Poisson and negative-binomial models require nonnegative integer counts. The count-likelihood runners report negative observations in calibration rows. [Observation mapping](forecasting_odemodels%20code/quantdiffObservationCurve.m)

## Rolling windows and forecast horizons

`tstart1` and `tend1` are the first and last **window-start row indices**, inclusive. To forecast from three successive origins, run the following from the code folder.

```matlab
opts = @options_forecast_SEIR_flu1918_IC;
rng(1, 'twister');
Run_Forecasting_ODEModel(opts, 1, 3, 17, 10);
plotForecast_ODEModel(opts, 1, 3, 17, 10);
```

The last window fits rows 3–19 and forecasts rows 20–29. With `getperformance = 1`, evaluation requires data through row `tend1 + windowsize1 - 1 + forecastingperiod`. If the runner shortens `tend1` to fit the data, use the reported end index in the plotter. For a projection beyond the observed series, set `getperformance = 0` in a copy of the forecast options. The saved `timevect2` grid supplies future time labels where CSV entries are `NaN`.

### Initial conditions for each window

In the quick start, the first observation sets the fitted state $C(0)$ and `vars.initial` supplies the other states. A later window using those options would restart with the same four infectious individuals as row 1. The runners issue `QuantDiffForecast:rollingInitialStates` when a window starts after row 1, the model has unfitted states and the options omit an initial-condition function.

The rolling example above uses `vars.initialFcn = @(p, x0, vars) ...` to initialize each window. Here `p` contains the current model parameters and `x0` contains the starting values of the fitted states. Return one finite, real, nonnegative value per state in `vars.label` order. The toolbox then sets the states indexed by `vars.fit_index` to `x0`. With `params.fixI0 = 1`, these entries equal the first observations in the window. With `params.fixI0 = 0`, they are estimated.

For the 1918 example, $I(0)$ is estimated, $E(0)$ is approximated by the first daily count divided by κ, and $R(0)$ is approximated by the number of cases reported before the window. The susceptible count is the remainder of the population. For your own model, `vars.window` provides the start row (`tstart`), calibration rows (`rows`), first time label (`time0`), first observations (`firstObs`) and the sum of each observed series over preceding rows (`priorCum`). The interpretation of that sum as removals is specific to this example. The example shows how to initialize later windows. It is not consistently the better forecaster for this series, as the comparison below shows. [1918 options](forecasting_odemodels%20code/options_forecast_SEIR_flu1918_IC.m) · [Initial-condition function](forecasting_odemodels%20code/quantdiffInitialConditions.m)

### Read a forecast score

With `getperformance = 1` and a fully observed forecast period, `Run_Forecasting_ODEModel` writes two score tables per fitted series and origin. WIS is the weighted interval score. Lower WIS indicates a better probabilistic forecast on the same targets.

| Evaluation question | File | Interpretation of row 10 |
| --- | --- | --- |
| How accurate was the forecast **on day 10**? | `scores-byhorizon-*.csv` | The observed count, predictive median, absolute and squared error, and WIS for day 10. `cover50` and `cover95` are 0/1 indicators for that observation. |
| How accurate were the forecasts **over days 1–10**? | `performance-forecasting-*.csv` | Mean absolute error (MAE), mean squared error (MSE), 95% prediction interval coverage (%) and mean WIS across the ten targets. The baseline comparison figure uses this cumulative WIS. |

The per-horizon table also reports `WIS_log`, computed after clipping negative values to zero and transforming observations and draws with `log1p`. When the runner can generate a flat baseline, `WIS_baseline` gives its score and `relWIS` gives the model-to-baseline ratio. For a positive baseline WIS, a ratio below 1 favors the model. [Per-horizon scoring](forecasting_odemodels%20code/quantdiffForecastScores.m) · [Cumulative point scores](forecasting_odemodels%20code/computeforecastperformance.m) · [Cumulative WIS](forecasting_odemodels%20code/computeWIS.m)

### Compare models across forecast origins

Compare the two SEIR initializations on the same observations. The rows of `windows` are `[first calibration row, number of calibration rows]`. These three expanding windows all begin at row 1. Each produces a ten-day forecast.

```matlab
models = {@options_forecast_SEIR_flu1918_dist1_3, @options_forecast_SEIR_flu1918_IC};
labels = {'Quick-start SEIR', 'Window-specific SEIR'};
windows = [1 14; 1 17; 1 20];
R = quantdiffCompareForecasts(models, windows, 10, 'Labels', labels, ...
    'Baseline', 'flat', 'Run', true, 'Seed', 1, 'Name', 'flu1918');
disp(R.byHorizon(R.byHorizon.horizon == 10, :))   % day 10 across origins
disp(R.summary)   % all shared targets
```

Read this comparison as a demonstration of the workflow. On this series neither initialization wins at every origin. The window-specific model fits the calibration rows more closely, yet from the 14- and 17-observation origins it forecasts worse than the quick-start model. From the 17-observation origin, its mean WIS over the ten held-out days is about 280, against 49 for the quick-start model, and its 95% prediction intervals contain none of the ten counts, against seven. From the 20-observation origin the ranking reverses, with a mean WIS of about 150 against 320. Averaged over the three origins, its mean WIS is about 181 against 132. A closer fit to the calibration rows does not imply a better forecast, and one origin is not enough to rank models.

The example reuses saved forecasts and runs any missing single-window forecasts. By default, `'Run', false` requires saved files. Give each model configuration a distinct `model.name` when changing settings, such as the bootstrap count, that do not appear in filenames. For files from a rolling run such as `Run_Forecasting_ODEModel(opts, 1, 3, 17, 10)`, use `[1 17; 2 17; 3 17]` and `'Tend', 3` to match the saved filenames. For multiple observed series, `'FitIndex', 2` selects the second fitted series.

`R.perHorizon` has one row per model, origin and target. `R.byHorizon` averages over origins for each lead, and `R.summary` averages over all targets with finite scores for every model and the chosen baseline. Summary `relWIS` is the mean model WIS divided by the mean baseline WIS. For daily counts, use `'Aggregate', 7` to score weekly totals. The scorer sums each predictive trajectory and the observations over successive seven-day blocks before scoring. A 14-day horizon gives two weekly targets. With a ten-day horizon, only the first seven days form a complete block. [Model comparison](forecasting_odemodels%20code/quantdiffCompareForecasts.m)

### Baselines and forecast ensembles

The flat baseline centers forecasts on the last calibration observation and samples historical changes. The trend baseline extrapolates a log-linear fit to the last seven observations by default. Select `'Baseline', 'flat'`, `'trend'` or `'none'` in the comparison above. For a single saved forecast, `quantdiffForecastScores` accepts the predictive draws (`forecast2`), the matching two-column observations (`[time, value]`), calibration length and horizon. Its optional `'BaselineDraws'` argument accepts draws from `baselineForecast`. [Baseline implementation](forecasting_odemodels%20code/baselineForecast.m)

<details>
<summary>Score a baseline directly for the quick-start forecast</summary>

```matlab
data = load(fullfile('input', 'curve-flu1918SF.txt'));
calib = data(1:17, :);
baselineDraws = baselineForecast(calib, 10, 'flat', 6000, 1);
[~, baselineWIS] = computeWIS(calib, data, baselineDraws, 10);
relativeWIS = scores.WIS(end) / baselineWIS(end, 2)
```

This uses the quick-start `scores` table and the same ten held-out counts. Substitute `'trend'` to compare with the trend baseline.

</details>

To combine the forecasts from the comparison above, reuse its saved 17-observation fits. This example gives the two models equal weight.

```matlab
E = run_ensemble(models, labels, -1, 17, 10, 'flu1918-equal', 'Plot', false);
disp(E.weights)
disp(E.relWIS_ens)   % relative WIS against the flat and trend baselines
```

The bundled pooling function combines member predictive distributions for the same series, window and horizon. It computes ensemble scores from the exact quantiles of the weighted empirical distribution. Weight codes are `-1` for equal weights, `0` or `1` for Akaike weights, `2` for inverse calibration WIS and `3` for softmax calibration WIS. [Ensemble runner](forecasting_odemodels%20code/run_ensemble.m) · [Pooling function](forecasting_odemodels%20code/getensemble_linearpool_from_curves.m) · [Weight definitions](forecasting_odemodels%20code/ensembleWeights.m)

## Profile likelihood

To ask how tightly the observed data constrain a parameter, fix it at successive values and reoptimize the remaining estimated quantities. This example profiles β over the first 17 observations and reestimates negative-binomial dispersion at each point. The profiler performs its own reference fit. [Profile-likelihood runner](forecasting_odemodels%20code/Run_ProfileLikelihood_ODEModel.m) · [Raue et al. (2009)](https://doi.org/10.1093/bioinformatics/btp358)

```matlab
rng(1, 'twister');
opts = struct('paramIndex', 1, 'alpha', 0.05, 'numstartpoints', 3);
PL = Run_ProfileLikelihood_ODEModel(@options_fit_SEIR_flu1918_dist1_3, 1, 17, opts);
disp(PL.referenceStatus)
disp(PL.param(1).status)
disp(PL.param(1).CI)
```

`opts.alpha` sets the significance level for the likelihood-ratio threshold. Read the approximate 95% confidence interval with the reference status and the profile plot. An endpoint at a search bound can reflect the chosen bound. If the optimizer failed or exhausted its search budget, inspect `stopReason` and the explored range before reporting an endpoint. The flatness index (`FI`) is the largest deviance at points free of flagged nuisance-parameter bounds divided by the threshold. A value near zero means the likelihood changed little across those points. Structural identifiability requires a separate analysis.

For a comparison across observation windows and choices of free parameters, run `plot_PL_1918_examples`. The script clears the workspace and produces a 12-panel figure from six configurations, including a simulated series with known parameters. In case C, both β and population size $N$ are estimated, and $S(0)$ follows each trial value of $N$. The first 17 observations favor a small population that is inconsistent with subsequent case counts. Here, interpreting the early slowdown as susceptible depletion gives a misleading estimate of population size. [Case C options](forecasting_odemodels%20code/options_fit_SEIR_flu1918_PL_C.m) · [Profile examples](forecasting_odemodels%20code/plot_PL_1918_examples.m)

## Practical identifiability simulations

To evaluate recovery under a specified observation design, simulate series with known parameter values, refit each series and check how often its bootstrap interval contains the truth. The bundled exponential-growth study estimates the growth rate `r` and dispersion α. Its generating values are `r = 0.18` and `α = 5`, with negative-binomial variance $\mu + 5\mu$. [Study script](forecasting_odemodels%20code/plotPracticalIdentifiabilityResults.m) · [Study options](forecasting_odemodels%20code/options_forecast_PII_EXPO_r_dist1_3.m)

Edit **USER SETTINGS** in `plotPracticalIdentifiabilityResults.m`, then run the script by name. It clears the workspace at startup. The shipped settings compare windows of 30, 40 and 50 observations, with 20 simulated datasets per window and 300 bootstrap refits per dataset. The replicate number seeds the simulation across window lengths, so comparisons by replicate are paired. For a short check of this workflow, `qdfSmokeTest3(1)` uses two replicates at one window length.

```matlab
plotPracticalIdentifiabilityResults
```

The replicate runner scores point errors, interval coverage and WIS from the predictive draws that include observation noise. Open `output/EXP_Negbin3-alpha5_W30-10-50/recovery_summary.csv` for the parameter recovery results. Each row represents one parameter and one window length.

| Question | Columns to read | Interpretation |
| --- | --- | --- |
| How many replicates contributed? | `files`, `converged` | Result files found and replicates eligible for summary. A documented failed point fit or mismatched noise record excludes a replicate. |
| Are the estimates close to the truth? | `truth`, `bias_MLE`, `RMSE_MLE`, `bias_bootmedian` | Bias measures systematic error. RMSE also reflects variability. Compare the fitted point estimate with the bootstrap median for evidence of skew or boundary effects. |
| Do the intervals cover the truth? | `coverage95`, `coverage_MCSE` | Empirical coverage is a **fraction**. With 20 eligible replicates, one missed interval changes coverage by 0.05. |
| How wide are the intervals? | `PII_scale`, `median_PII`, `median_rangeFrac`, `frac_bound_hits` | Interpret width on the selected scale and check how often estimates accumulate at optimizer bounds. |

In the default study, **read `median_rangeFrac` for `r`** because its optimizer bounds are [0, 10]. For α, the lower bound is positive, so read `median_PII`. Use `run_flag = false` to summarize existing files. Changing truth, constraints or bootstrap settings calls for a distinct results folder and simulated-input base name (`cadfilename1`), since the script skips completed window/replicate pairs. The input files have companion `-noise.mat` records. For Poisson and likelihood-based negative-binomial studies, the summary checks available records against the requested settings. Older results with missing records can enter the summary with a warning. [Recovery summary implementation](forecasting_odemodels%20code/summarizeIdentifiability.m)

<details>
<summary>Interval-width measures and simulation options</summary>

For bootstrap interval endpoints $L,U$ and optimizer bounds $L_b,U_b$, the reported measures are as follows.

| Measure | Definition |
| --- | --- |
| Parameter identifiability index (PII), stored as `SCI` | $\log_{10}(U/L)$ when both endpoints are positive. A value of 1 means the interval spans a factor of ten. `n_PII_undefined` counts replicates without a defined log width. |
| `PII_rangeFrac` | $(U-L)/(U_b-L_b)$, clipped to [0, 1] for a free parameter. A value of 0.8 means the interval covers 80% of its allowed linear range. |
| `PII_box`, `PII_ratio` | The log span of positive optimizer bounds and the interval’s share of that span. On the range-fraction scale, `PII_ratio` equals `PII_rangeFrac`. |
| `fracLB`, `fracUB` | Fractions of bootstrap estimates at the lower and upper bounds, using a tolerance of 0.1% of each bound’s magnitude. `frac_bound_hits` counts the fraction of replicates with `fracLB + fracUB > 0.05`. |

Set `params.PIIscale` by parameter slot to `1` (log width), `0` (range fraction) or `NaN` (automatic). Read width alongside coverage and bound hits. Bootstrap estimates concentrated at one bound can give a narrow interval with poor coverage. [Index definitions](forecasting_odemodels%20code/quantdiffPII.m) · [Bound-hit definitions](forecasting_odemodels%20code/boundHitFractions.m)

To simulate one dataset without fitting, or to fit selected replicates directly, use a forecast-type options function.

```matlab
studyOpts = @options_forecast_PII_EXPO_r_dist1_3;
rng(1, 'twister');
[states, observations, noiseTruth] = Run_simulate_ODEModel(studyOpts, 30, 5);
outDir = fullfile(pwd, 'output', 'EXP_Negbin3-alpha5_W30');
Run_PracticalIdentifiability_ODEModel(studyOpts, 30, 5, 2, [], outDir);
```

The simulator returns noise-free states, a matrix of time and noisy observations, and the generating noise settings. For method 5, choose the bundled method-5 options and pass `[alpha d]` as the noise factor for variance $\mu+\alpha\mu^d$. A scalar sets $d=1$. An empty noise factor uses `[0.5 1.5]` for method 5. [Simulator](forecasting_odemodels%20code/Run_simulate_ODEModel.m) · [Replicate runner](forecasting_odemodels%20code/Run_PracticalIdentifiability_ODEModel.m) · [Method-5 options](forecasting_odemodels%20code/options_forecast_PII_EXPO_r_dist1_5.m)

</details>

## Configuration

The fitting and profile runners expect a fit-type options function (13 outputs). Forecasting, simulation and parameter recovery use a forecast-type function (15 outputs). Keep the output order when you copy a template. Explicit window and horizon arguments in a runner call override the corresponding options.

| Setting | What to specify |
| --- | --- |
| `cadfilename1`, `caddisease`, `datatype` | Input filename, dataset label and observed quantity. The labels appear in output filenames. |
| `model.fc`, `model.name` | ODE function handle and a label for this model or analysis. |
| `params.label`, `params.initial`, `params.LB`, `params.UB` | Parameter names, starting values and finite bounds in the order used by the ODE. Starting values must lie within the bounds. |
| `params.fixed`, `params.fixI0` | `fixed = 1` holds the corresponding parameter at its initial value. `fixI0 = 1` pins each fitted state’s initial value to its first observation. Set a parameter’s `fixed` entry to `0` to estimate it, or set `fixI0 = 0` to estimate the fitted states’ initial values. Bootstrap datasets retain the first observation of each series in both cases. |
| `vars.label`, `vars.initial`, `vars.fit_index`, `vars.fit_diff` | State names, initial conditions, observed-state indices and level/increment flags. |
| `vars.initialFcn` (optional) | Function of the current parameters, fitted initial values and calibration window. See [Initial conditions for each window](#initial-conditions-for-each-window). |
| `params.composite`, `params.composite_name`, `params.extra0` | Optional derived-quantity function and label, and extra data passed to the ODE. Set the function to `[]` when unused. |
| `method1`, `dist1` | Fitting objective and observation-noise model. The options function sets the global `method1`. See the combinations below. |
| `B`, `numstartpoints` | Bootstrap datasets per window and starting-point budget for the original fit. The optimizer also adds seed and perturbed starts. |
| `windowsize1`, `forecastingperiod`, `tstart1`, `tend1` | Calibration length, forecast length, first window-start row and last window-start row. |
| `getperformance`, `printscreen1` | Enable forecast evaluation/export and selected displays, respectively. Some scoring and console output also occur with these flags off. |

### Estimation methods and observation models

| `method1` | `dist1` | Fitting objective | Noise for bootstrap and prediction |
| ---: | ---: | --- | --- |
| 0 | 0 | Least squares | Gaussian |
| 0 | 1 | Least squares | Poisson |
| 0 | 2 | Least squares | Negative binomial with variance `factor1 * mu` estimated from Pearson residuals |
| 1 | 1 | Poisson likelihood | Poisson with variance $\mu$ |
| 3 | 3 | Negative-binomial likelihood | Variance $\mu+\alpha\mu$ |
| 4 | 4 | Negative-binomial likelihood | Variance $\mu+\alpha\mu^2$ |
| 5 | 5 | Negative-binomial likelihood | Variance $\mu+\alpha\mu^d$ |
| 6 | 6 | Sum of absolute deviations | Laplace with scale estimated from absolute residuals |

Here $\mu$ is the predicted observation, $\alpha$ is dispersion and $d$ is the variance exponent. For positive method codes, the runners set `dist1 = method1`. For method 0, all three noise choices produce the same unweighted least-squares objective. Code 2 is available only for `dist1`. If its estimated variance-to-mean factor is unavailable or at most one, that window uses Poisson draws and retains `dist-2` in filenames. Gaussian and Laplace noise can yield negative draws. The runners warn when the expected negative fraction exceeds 5% and record it in `noiseLog`. [Pearson dispersion](forecasting_odemodels%20code/estimatePearsonDispersion.m) · [Laplace scale](forecasting_odemodels%20code/estimateLaplaceScale.m) · [Negative-draw diagnostic](forecasting_odemodels%20code/quantdiffNegativeNoiseShare.m)

<details>
<summary>Optimizer settings and AICc</summary>

The named quick-start options and the default forecast options enable parallel local fits, log-scale search and polishing. Set `params.useParallel = false` in a copy for serial fitting. If a parallel pool is unavailable, the optimizer falls back to serial execution. With `params.logscale = true`, free coordinates with positive lower bounds are searched on a log scale and reported in their original units.

| Optional field | Default when omitted | Purpose |
| --- | --- | --- |
| `params.useParallel`, `params.logscale` | `false`, `false` | Parallel local fits and log-scale search, respectively |
| `params.optimOptions` | `[]` | Override selected `fmincon` settings |
| `params.polish` | `false` | Try Nelder–Mead refinement and retain an improved objective |
| `params.polish_gradient_tol`, `params.polish_budget` | `1e-2`, `100` | Refinement gradient threshold and evaluation budget per free parameter |

For $s$ observed series and $W$ rows per series, the effective sample size is $s(W-1)$ when initial observations are pinned and $sW$ when the initial values are estimated. The corrected Akaike information criterion (AICc) requires at least $k+2$ informative observations for $k$ estimated quantities. The quick start has 16 informative observations and two estimated quantities (β and α). Least-squares AICc includes one Gaussian variance regardless of bootstrap noise. Method 6 uses the Laplace likelihood and counts its scale. Count-likelihood AICc restores the data-dependent log-factorial term and removes the pinned first observations’ contribution. Compare AICc values for the same responses, calibration rows and initial-value treatment. [Optimizer](forecasting_odemodels%20code/fit_model.m) · [AICc](forecasting_odemodels%20code/getAICc.m) · [Effective sample size](forecasting_odemodels%20code/quantdiffEffectiveN.m)

</details>

## Outputs

All paths below are relative to `forecasting_odemodels code/output/`. A forecast trajectory is saved for each observed series and window. For parameter estimates, start with the parameter CSV and then check the corresponding MAT file for convergence and bound diagnostics.

| File or prefix | What it contains |
| --- | --- |
| `Forecast-model_name-*.csv` | Observations, predictive median, 95% prediction limits and the runner’s `model_median`. The forecast plotter rewrites the file with five columns and uses the raw predictive median. The runner clips its exported medians and limits at zero. |
| `performance-calibration-*.csv`, `performance-forecasting-*.csv` | MAE, MSE, 95% prediction interval coverage (%) and WIS. Row `h` of the forecasting file summarizes leads 1 through `h` at one origin. |
| `scores-byhorizon-*.csv` | One forecast lead per row for a single model and origin. Includes individual errors, coverage indicators, WIS, log-scale WIS and, when available, flat-baseline scores. |
| `scores-comparison-<name>-perhorizon.csv`, `-byhorizon.csv`, `-summary.csv` | Model comparison at each target, by lead across origins, and across all shared targets. `quantdiffCompareForecasts` writes these by default. |
| `quantile-*.csv`, `StateVars-*.csv` | The plotter exports 23 predictive quantiles from 0.01 to 0.99. State files summarize bootstrap trajectories. In rolling runs, successive windows reuse the same state CSV filename. Run one start at a time to keep one state CSV per origin. |
| `parameters-rollingwindow-*.csv`, `parameters-composite-*.csv` | Point estimates from the observed data and 2.5th/97.5th bootstrap percentiles. For a derived quantity, the central estimate is its value at the observed-data fit. |
| `AICc-*.csv` | AICc, its two additive components and the count of estimated quantities. |
| `SCIs-rollingwindow-*.csv`, `PII*-rollingwindow-*.csv` | Log width, optimizer-box width, normalized width, range fraction, selected scale and bound-hit fractions. |
| `MCSEs-rollingwindow-*.csv` | Monte Carlo standard errors of interval endpoints and log width, estimated by resampling the bootstrap estimates 500 times. `_SE_bootmean` is the standard error of the bootstrap mean. |
| `bootstraps-ODEModel-*.mat`, `parameters-ODEModel-*.mat` | Bootstrap estimates and noise details in the former. Parameter summaries, `param_bootmedian`, `estimOutsideCI` and optimizer diagnostics (`numericalAudit`) in the latter. |
| `Forecast-ODEModel-*.mat` | Saved forecast for a window and fitted series. `forecast1` contains bootstrap model trajectories for that series. `forecast2` includes observation noise. |
| `PLsummary-*.csv`, `PLprofile-param-*.csv`, `PL-*.mat` | Profile intervals and status, pointwise profiles and the full profile result. The default run also saves a PNG. |
| `results_W<W>_rep<r>.mat`, `recovery_summary.csv` | Results for each simulated replicate and the recovery summary in the study folder. |

The runners report requested, accepted and unresolved bootstrap fits. Unresolved refits can still enter interval estimates, so inspect `numericalAudit` before interpreting narrow intervals. Save the seed, options file, source revision and outputs together. Changing `B`, parameter bounds or optimizer settings can reuse the same filenames.

## Adding a model

Place the ODE function in the code folder. It takes four inputs (time, states, parameters and `params.extra0`) and returns a column vector of derivatives in state order. The solver enforces nonnegative states. The bundled exponential-growth model is a complete one-state example. [EXP.m](forecasting_odemodels%20code/EXP.m) · [SEIR1.m](forecasting_odemodels%20code/SEIR1.m) · [Solver settings](forecasting_odemodels%20code/quantdiffODEOptions.m)

For example, copy a quick-start options function, rename it, and replace the model-specific assignments with the following to fit increments of $C$ under $dC/dt=rC$ using the existing `EXP.m`. Retain the data, estimation and window settings from the template. For your own equations, replace `@EXP` with a handle to your four-input function.

```matlab
model.fc = @EXP;
model.name = 'Exponential growth';
params.label = {'r'};
params.initial = 0.1;
params.LB = 0;
params.UB = 1;
params.fixed = 0;
params.fixI0 = 1;
params.composite = [];
params.composite_name = '';
params.extra0 = [];

vars.label = {'C'};
vars.initial = 4;
vars.fit_index = 1;
vars.fit_diff = 1;
```

A composite function receives a matrix with one row per parameter draw and returns one value per row. The bundled SEIR example computes $R_0=\beta/\gamma$ across those rows. [Composite example](forecasting_odemodels%20code/R0s.m)

## Example visualizations

The first four figures show an example fit to the 1918 San Francisco influenza series. In the fit and forecast panels, the first 17 observations cover days 0–16 and the ten held-out observations cover days 17–26. In plots labeled $C'(t)$, a plotted value represents the daily increment of the cumulative case state. These screenshots illustrate the toolbox displays. The four comparisons below have a separate regeneration script.

<table>
  <tr>
    <td width="50%" align="center">
      <img src="docs/images/model_fit.png" alt="Example SEIR fit to the observed time series" width="100%"><br>
      <sub>Fit to daily counts on days 0–16. Blue circles are observations. Red dashed lines enclose the central 95% prediction interval.</sub>
    </td>
    <td width="50%" align="center">
      <img src="docs/images/forecast.png" alt="Example model forecast with uncertainty intervals" width="100%"><br>
      <sub>Ten-day forecast. The vertical dashed line marks the last calibration day (day 16).</sub>
    </td>
  </tr>
  <tr>
    <td width="50%" align="center">
      <img src="docs/images/parameters.png" alt="Example estimated parameters and interval summaries" width="100%"><br>
      <sub>Bootstrap distribution of β. The red dashed line marks the fit to the observed counts. κ, γ and N are fixed.</sub>
    </td>
    <td width="50%" align="center">
      <img src="docs/images/R0.png" alt="Example bootstrap distribution of the basic reproduction number" width="100%"><br>
      <sub>Bootstrap distribution of the derived reproduction number, R₀ = β/γ, with γ fixed at 1/4.1 per day.</sub>
    </td>
  </tr>
</table>

The next four figures use the shipped 1918 options and the same input series. From the code folder, run `make_readme_figures(true)` to recompute five forecasts and one profile likelihood and write these four PNGs to `docs/images/`. Calling `make_readme_figures` without an argument reuses completed results in `output/readme_figures_cache.mat`. [Figure script](forecasting_odemodels%20code/make_readme_figures.m)

<p align="center">
  <img src="docs/images/forecast_by_calibration_length.png" alt="Ten-day SEIR forecasts from the first 14, 17 and 20 observations with 50% and 95% prediction intervals" width="100%"><br>
  <sub>Ten-day forecasts after 14, 17 and 20 observations. Each window has a different set of ten held-out days. Shading shows the 50% and 95% prediction intervals (<a href="#rolling-windows-and-forecast-horizons">forecast evaluation</a>)</sub>
</p>

<p align="center">
  <img src="docs/images/noise_model_comparison.png" alt="The same SEIR forecast under least squares with normal noise, Poisson likelihood and negative-binomial likelihood" width="100%"><br>
  <sub>Three methods fitted to days 0–16 and scored on days 17–26. The 95% prediction intervals contain 2 of 10 counts for least squares with Gaussian noise, 6 of 10 for Poisson, and 7 of 10 for negative binomial (<a href="#estimation-methods-and-observation-models">methods</a>)</sub>
</p>

<table>
  <tr>
    <td width="62%" align="center">
      <img src="docs/images/baseline_comparison.png" alt="SEIR forecast compared with the flat and trend baselines, with weighted interval score by forecast horizon" width="100%"><br>
      <sub>At lead h, each curve shows mean WIS over leads 1 through h. Over all ten days the SEIR-to-baseline WIS ratios are 0.30 for trend and 0.11 for flat (<a href="#baselines-and-forecast-ensembles">baselines</a>)</sub>
    </td>
    <td width="38%" align="center">
      <img src="docs/images/profile_likelihood_beta.png" alt="Profile likelihood of the transmission rate with the 95% threshold and interval" width="100%"><br>
      <sub>Profile likelihood of β using the first 17 observations. The shaded area is the approximate 95% likelihood-ratio interval for β (<a href="#profile-likelihood">profile likelihood</a>)</sub>
    </td>
  </tr>
</table>

<details>
<summary>Additional model-state and forecast-performance illustrations</summary>

**Before fitting.** The blue curves use sampled parameter values. The green curve has the smallest sum of squared errors among those draws. The right panel zooms in on the observations. [Model exploration plotter](forecasting_odemodels%20code/plotODEModel.m)

![Example model trajectories before fitting](docs/images/model_solutions.png)

**State trajectories.** The five panels show S, E, I, R and cumulative cases C. The solid line is the bootstrap median and the dashed lines are the 2.5th and 97.5th percentiles of the model states. [Fit plotter](forecasting_odemodels%20code/plotFit_ODEModel.m)

![Example state-variable trajectories and uncertainty summaries](docs/images/stateVars.png)

**Forecast performance from an earlier run.** Each point summarizes leads 1 through h. The 95% coverage curve ends at 90% (9 of 10 days). The scripted 17-observation negative-binomial comparison above shows 70% (7 of 10). Use the score CSV from the same run as the forecast when reporting coverage or WIS. [Forecast plotter](forecasting_odemodels%20code/plotForecast_ODEModel.m)

![Cumulative MAE, MSE, 95% prediction interval coverage and WIS from an earlier forecast](docs/images/forecastingPerformance.png)

</details>

## Checks and support

| Task or symptom | Command or remedy |
| --- | --- |
| Check input validation | `qdfSmokeTest([7 8])` checks that unsupported methods and negative count observations are caught. |
| Exercise the main workflows | `qdfSmokeTest()` runs the first eight tests. `qdfSmokeTest2()` adds ten checks, including multiple series and rolling windows. `qdfSmokeTest3(1:3)` covers recovery plotting, a fit with three free model parameters and a profile of κ. `qdfSmokeTest3(4)` fits and plots the bundled Spain under-reporting example using the first 60 observations. |
| Check initial conditions without fitting | `qdfTestInitialConditions()` checks state initialization and objective values against reference values. [Initial-condition checks](forecasting_odemodels%20code/qdfTestInitialConditions.m) |
| Check forecast scores without fitting | `qdfTestForecastScores()` checks per-lead and cumulative scores, block totals, baseline ratios and missing targets. [Score checks](forecasting_odemodels%20code/qdfTestForecastScores.m) |
| MATLAB cannot find a function | Run `which Run_Forecasting_ODEModel -all` and `ver` to check the active source folder and installed products. |
| The plotter cannot find a forecast | Use the options and effective window indices, calibration length and horizon from the completed forecast. |
| The options function returns too few outputs | Use a fit-type template for fitting or profiling and a forecast-type template for forecasting or simulation. |

The smoke tests print a verdict for each requested workflow. [First test suite](forecasting_odemodels%20code/qdfSmokeTest.m) · [Second suite](forecasting_odemodels%20code/qdfSmokeTest2.m) · [Third suite](forecasting_odemodels%20code/qdfSmokeTest3.m)

For a software issue, include the command, options file, source revision, MATLAB version, error text and a small shareable dataset. For scientific questions, contact Gerardo Chowell at Georgia State University. [GitHub Issues](https://github.com/gchowell/QuantDiffForecast-Toolbox/issues) · [gchowell@gsu.edu](mailto:gchowell@gsu.edu)

## Citation

Chowell G, Bleichrodt A, Luo R. Parameter estimation and forecasting with quantified uncertainty for ordinary differential equation models using QuantDiffForecast: A MATLAB toolbox and tutorial. *Statistics in Medicine*. 2024, 43(9), 1826–1848. [DOI](https://doi.org/10.1002/sim.10036) · [PubMed](https://pubmed.ncbi.nlm.nih.gov/38378161/)

<details>
<summary>BibTeX</summary>

```bibtex
@article{Chowell2024QuantDiffForecast,
  author  = {Gerardo Chowell and Amanda Bleichrodt and Ruiyan Luo},
  title   = {Parameter estimation and forecasting with quantified uncertainty
             for ordinary differential equation models using {QuantDiffForecast}:
             A {MATLAB} toolbox and tutorial},
  journal = {Statistics in Medicine},
  year    = {2024},
  volume  = {43},
  number  = {9},
  pages   = {1826--1848},
  doi     = {10.1002/sim.10036}
}
```

</details>

For background on fitting dynamic models and quantifying uncertainty, see [Chowell (2017)](https://doi.org/10.1016/j.idm.2017.08.001).

## License

This toolbox is distributed under CC0 1.0 Universal. [License](LICENSE)
