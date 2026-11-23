% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

function [cadfilename1,caddisease,datatype, dist1, numstartpoints,B, model, params,vars,getperformance, forecastingperiod,windowsize1,tstart1,tend1,printscreen1]=options_forecast_SEIR_flu1918_IC

% OPTIONS_FORECAST_SEIR_FLU1918_IC  SEIR forecast with I(0) estimated and initial conditions that follow the window
%
% Example of vars.initialFcn (see quantdiffInitialConditions). It differs from
% options_forecast_SEIR_flu1918_dist1_3 in three ways:
%   1. I_0 is parameter 5 and is estimated with beta. SEIR1 reads parameters
%      1-4 only, so the model function is unchanged.
%   2. vars.initialFcn builds the initial conditions from the parameters and
%      the data: E(0) = (first observation of the window)/kappa, because new
%      cases per day are about kappa*E, and
%      S(0) = N - I_0 - E(0) - (cases reported before the window).
%   3. A calibration window that starts after row 1 therefore starts from a
%      depleted population and its own I(0) and E(0), not from vars.initial.
%   E_0 (parameter 6) is not read by this function and stays fixed. An
%   alternative vars.initialFcn that also estimates E_0 is given below; with 17
%   observations I_0 and E_0 are then only weakly identified separately (their
%   bootstrap intervals reach the lower bound).
%
% Overview
%   Returns all configuration needed by Run_Forecasting_ODEModel to (i) calibrate
%   an ODE model over rolling windows and (ii) generate out-of-sample forecasts
%   with quantified uncertainty (parametric bootstrap).
%
% Usage
%   [cadfilename1, caddisease, datatype, dist1, numstartpoints, B, ...
%    model, params, vars, getperformance, forecastingperiod, ...
%    windowsize1, tstart1, tend1, printscreen1] = options_forecast_SEIR_flu1918_IC;
%
% Typical runs
%   Run_Forecasting_ODEModel(@options_forecast_SEIR_flu1918_IC, 1, 1, 17, 10);     % first window
%   Run_Forecasting_ODEModel(@options_forecast_SEIR_flu1918_IC, 21, 21, 17, 10);   % window that starts at row 21
%
% Outputs (returned by this options function)
%   cadfilename1     char     Base filename of input series (looked for in ./input)
%   caddisease       char     Label used in figures/outputs (e.g., '1918 Flu')
%   datatype         char     'cases' | 'deaths' | 'hospitalizations' | ...
%   dist1            scalar   Error-structure code (see “Estimation & error models”)
%   numstartpoints   int      MultiStart initial guesses per optimization
%   B                int      Bootstrap replicates for CI/PI and forecast uncertainty
%   model            struct   Model definition:
%                              .fc    (function handle) RHS ODE, e.g., @SEIR1
%                              .name  (char)           Human-readable model name
%   params           struct   Parameter settings:
%                              .num       (int)      number of parameters
%                              .label     (cellstr)  symbols (e.g., {'\beta','\kappa','\gamma','N'})
%                              .LB,.UB    (1×num)    lower/upper bounds
%                              .initial   (1×num)    starting values
%                              .fixed     (1×num)    1=fixed at .initial, 0=estimated
%                              .fixI0     (logical)  1=fix initial observed state to first datum; 0=estimate
%                              .composite (function) handle for composite metric (e.g., @R0s) or []
%                              .composite_name (char) name of composite metric (e.g., 'R_0')
%                              .extra0    (any)      optional constants passed into the model function
%                             Optimizer settings (read by fit_model; all optional):
%                              .logscale  (logical)  1=search in log10 of positive-bounded parameters
%                                                    (recommended when an estimated parameter spans
%                                                    several orders of magnitude, e.g. K, N, alpha)
%                              .optimOptions (struct) fmincon overrides, e.g.
%                                                    struct('StepTolerance',1e-8,'FunctionTolerance',1e-8);
%                                                    [] keeps the toolbox defaults
%                              .useParallel (logical) 1=run MultiStart local solves on a parallel pool
%                                                    (Parallel Computing Toolbox; serial if no pool)
%                              .polish    (logical)  1=Nelder-Mead polish of the MultiStart solution
%                              .polish_gradient_tol (scalar) polish only if the scaled gradient at the
%                                                    sqp solution exceeds this value (default 1e-2)
%                              .polish_budget (int)  Nelder-Mead evaluations per free parameter (default 100)
%   vars             struct   State/fit settings:
%                              .num      (int)       number of state variables
%                              .label    (cellstr)   names (e.g., {'S','E','I','R','C'})
%                              .initial  (1×num)     initial conditions
%                              .fit_index(int)       index of state matched to data (e.g., 5 for 'C')
%                              .fit_diff (logical)   1=fit derivative/incidence; 0=fit level
%   getperformance   logical  1=compute forecast performance metrics; 0=skip
%   forecastingperiod int     Forecast horizon (time steps ahead)
%   windowsize1      int      Rolling-window length (time steps)
%   tstart1          int      Start index of rolling-window analysis
%   tend1            int      End index of rolling-window analysis
%   printscreen1     logical  1=show figures/console progress; 0=silent
%
% Data format (./input/<cadfilename1>.txt, default)
%   Two columns with no header:
%     Col 1: time index (0,1,2,...)
%     Col 2: observed counts (nonnegative integers)
%
% Estimation & error models
%   A global variable 'method1' selects the estimation method; 'dist1' selects
%   the observational error and is synchronized to 'method1' as follows:
%
%     method1 = 0  Nonlinear least squares (LSQ). The objective is the same
%         unweighted sum of squares for every dist1; dist1 only selects the
%         noise used for the bootstrap and the prediction intervals:
%           dist1=0  Normal noise (constant variance)
%           dist1=1  Poisson noise (variance = mean)
%           dist1=2  Negative binomial noise with var = factor1·mean
%                    (factor1 estimated from the Pearson residuals of the fit)
%
%     method1 = 1  MLE with Poisson errors          => dist1=1 (set automatically)
%     method1 = 3  MLE NegBin: var = mean + α·mean  => dist1=3 (set automatically)
%     method1 = 4  MLE NegBin: var = mean + α·mean^2=> dist1=4 (set automatically)
%     method1 = 5  MLE NegBin: var = mean + α·mean^d=> dist1=5 (set automatically)
%     method1 = 6  Sum of absolute deviations       => dist1=6 (Laplace noise, set automatically)
%
% Notes
%   • MultiStart with a sensible 'numstartpoints' helps avoid local minima.
%   • Set params.fixI0=1 to anchor the initial observed state to the first datum;
%     set to 0 to estimate it jointly with other parameters.
%   • (tstart1, tend1) index rolling windows, not calendar dates.
%   • Defaults in this file instantiate an SEIR model for the 1918 SF dataset.
%   • Convergence of every fit can be checked in fitDiagnostics (returned by
%     fit_model): startAgreement close to 1 and a small scaledGradientInfNorm
%     mean the local runs agree and the solution is stationary.
%
% Dependencies
%   Optimization Toolbox (e.g., fmincon as local solver)
%   Global Optimization Toolbox (MultiStart)
%
% Author
%   Gerardo Chowell (Georgia State University) | gchowell@gsu.edu
%
% See also
%   Run_Forecasting_ODEModel, Run_Fit_ODEModel, plotForecast_ODEModel, plotFit_ODEModel



% <============================================================================>
% <=================== Declare global variables =======================================>
% <============================================================================>

global method1 % Parameter estimation method

% <============================================================================>
% <================================ Datasets properties =======================>
% <============================================================================>
% Located in the input folder, the time series data file is a text file with extension *.txt. 
% The time series data file contains the incidence curve of the epidemic of interest. 
% The first column corresponds to time index: 0,1,2, ... and the second
% column corresponds to the observed time series data.

cadfilename1='curve-flu1918SF'; % Name of the data file (without extension) containing the incidence curve

caddisease='1918 Flu'; % string indicating the name of the disease related to the time series data

datatype='cases'; % string indicating the nature of the data (cases, deaths, hospitalizations, etc)

% <=============================================================================>
% <=========================== Parameter estimation ============================>
% <=============================================================================>

method1=3; % Type of estimation method
% The negative binomial is kept because this file is also meant for later
% windows. On rows 1-17 the dispersion alpha goes to its lower bound, so
% method1=1 (Poisson) gives the same estimates with one estimated quantity
% fewer and a lower AICc; on rows 21-37 the data are strongly overdispersed
% and a Poisson fit is far worse.

% Nonlinear least squares (LSQ)=0,
% MLE Poisson=1,
% MLE (Neg Binomial)=3, with VAR=mean+alpha*mean;
% MLE (Neg Binomial)=4, with VAR=mean+alpha*mean^2;
% MLE (Neg Binomial)=5, with VAR=mean+alpha*mean^d;
% SAD (sum of absolute deviations) = 6.


dist1=3; % Define dist1 which is the type of error structure. See below:

%dist1=0; % Normal distribution to model error structure (method1=0)
%dist1=1; % Poisson error structure (method1=0 OR method1=1)
%dist1=2; % Neg. binomial error structure where var = factor1*mean where
                  % factor1 is empirically estimated from the time series
                  % data (method1=0)
%dist1=3; % MLE (Neg Binomial) with VAR=mean+alpha*mean  (method1=3)
%dist1=4; % MLE (Neg Binomial) with VAR=mean+alpha*mean^2 (method1=4)
%dist1=5; % MLE (Neg Binomial)with VAR=mean+alpha*mean^d (method1=5)
%dist1=6; % SAD method, Laplace distribution (method1=6).


switch method1
    case 1
        dist1=1;
    case 3
        dist1=3;
    case 4
        dist1=4;
    case 5
        dist1=5;
    case 6
        dist1=6;
end


numstartpoints=10; % Number of initial guesses for optimization procedure using MultiStart

B=300; % number of bootstrap realizations to characterize parameter uncertainty

% <==============================================================================>
% <============================== ODE model =====================================>
% <==============================================================================>
model.fc=@SEIR1; % name of the model function
model.name='SEIR model IC';   % string indicating the name of the ODE model (its own name, so its output files do not replace those of the other SEIR examples)

params.label={'\beta','\kappa','\gamma','N','I_0','E_0'};  % list of symbols to refer to the model parameters (I_0, E_0: initial infectious and exposed; read by vars.initialFcn only, and E_0 only by the alternative function below)
params.LB=[0.001 0 0 20 0.01 0.01]; % lower bound values of the parameter estimates
params.UB=[10 2 2 1000000 100000 100000]; % upper bound values of the parameter estimates
params.initial=[0.76 1/1.9 1/4.1 550000 4 4]; % initial parameter values/guesses

params.fixed=[0 1 1 1 0 1]; % Boolean vector to indicate any parameters that should remain fixed (1) to initial values indicated in params.initial. Otherwise the parameter is estimated (0). Here beta and I_0 are estimated; E_0 is not read by the function below.
params.fixI0=1; % Boolean variable indicating if the initial value of the fitting variable is fixed according to the first observation in the time series (1). Otherwise, it will be estimated along with other parameters (0).
params.composite=@R0s;  % Estimate a composite function of the individual model parameter estimates otherwise it is left empty.
params.composite_name='R_0'; % Name of the composite parameter
params.extra0=[]; % optional constants passed into the model function (empty if none)

% <==============================================================================>
% <========================== Optimizer settings (fit_model) ====================>
% <==============================================================================>
% All optional; see the header of fit_model.m. The optimizer searches only the
% estimated parameters, each rescaled to [0,1] by its bounds.

params.logscale=true; % search in log10 of positive-bounded parameters (recommended when an estimated parameter spans several orders of magnitude, e.g. K, N or alpha in methods 3-5)
params.optimOptions=[]; % fmincon overrides, e.g. struct('StepTolerance',1e-8,'FunctionTolerance',1e-8,'FiniteDifferenceStepSize',1e-5); [] keeps the toolbox defaults
params.useParallel=true; % run the MultiStart local solves on a parallel pool (Parallel Computing Toolbox); falls back to serial if no pool is available
params.polish=true; % Nelder-Mead polish of the MultiStart solution when the gradient test below says it is not stationary
params.polish_gradient_tol=1e-2; % run the polish only if the scaled gradient at the sqp solution exceeds this value
params.polish_budget=100; % Nelder-Mead evaluations per free parameter

vars.label={'S','E','I','R','C'}; % list of symbols to refer to the variables included in the model
vars.initial=[params.initial(4)-params.initial(5)-4/params.initial(2), 4/params.initial(2), params.initial(5), 0, 4];  % initial conditions at params.initial for a first observation of 4 (what the function below returns): the fallback, and the seed of the fitted state C when params.fixI0=0

% Initial conditions as a function of the parameters p, the initial value x0 of
% the fitted state C (the first observation of the window when params.fixI0=1)
% and the calibration window. States: S, E, I, R, C.
%   E(0) = x0/kappa: new cases per day are about kappa*E.
%   I(0) = I_0, estimated.
%   R(0) = vars.window.priorCum, the cases reported before the window (0 for a
%          window that starts at row 1), counted as removed: an approximation.
%   S(0) = N minus the three above.
vars.initialFcn=@(p, x0, vars) [p(4)-p(5)-x0/p(2)-vars.window.priorCum, x0/p(2), p(5), vars.window.priorCum, x0];

% Alternative that also estimates E(0), as parameter 6. To use it, make TWO
% changes: comment the line above and uncomment the line below, AND set
% params.fixed=[0 1 1 1 0 0]. (With the line below and params.fixed(6)=1,
% E(0) stays at params.initial(6) in every window.) I_0 and E_0 are then only
% weakly identified separately.
%vars.initialFcn=@(p, x0, vars) [p(4)-p(5)-p(6)-vars.window.priorCum, p(6), p(5), vars.window.priorCum, x0];

vars.fit_index=5; % index of the model's variable that will be fit to the observed time series data
vars.fit_diff=1; % boolean variable to indicate if the derivative of model's fitting variable should be fit to data.

% <==============================================================================>
% <========================== Forecasting parameters ===================================>
% <==============================================================================>

getperformance=1; % flag or indicator variable (1/0) to calculate forecasting performance or not

forecastingperiod=10; % forecast horizon (number of time units ahead)

% <==================================================================================>
% <========================== Parameters of the rolling window analysis =========================>
% <==================================================================================>

windowsize1=17;  % moving window size

tstart1=1; % time point for the start of rolling window analysis

tend1=1;  %time point for the end of the rolling window analysis

printscreen1=1;
