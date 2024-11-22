% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

function [cadfilename1, caddisease, datatype, dist1, numstartpoints, B, model, params, vars, windowsize1, tstart1, tend1, printscreen1] = options_fit_SEIR_flu1918

% <============================================================================>
% <=================== Declare Global Variables ==============================>
% <============================================================================>
% Declare global variables used in this function.
global method1; % Parameter estimation method

% <============================================================================>
% <========================= Dataset Properties ==============================>
% <============================================================================>
% The time series data file contains the incidence curve of the epidemic of interest.
% - The first column corresponds to the time index (e.g., 0, 1, 2, ...).
% - The second column contains the observed time-series data.

cadfilename1 = 'curve-flu1918SF'; % Name of the time-series data file
caddisease = '1918 Flu';          % Name of the disease related to the time-series data
datatype = 'cases';               % Nature of the data (e.g., cases, deaths, hospitalizations)

% <============================================================================>
% <======================= Parameter Estimation ==============================>
% <============================================================================>
% Define the method for parameter estimation and associated error structure.

method1 = 3; % Parameter estimation method:
% 0 - Nonlinear least squares (LSQ)
% 1 - Maximum Likelihood Estimation (MLE) Poisson
% 3 - MLE Negative Binomial (VAR = mean + alpha * mean)
% 4 - MLE Negative Binomial (VAR = mean + alpha * mean^2)
% 5 - MLE Negative Binomial (VAR = mean + alpha * mean^d)
% 6 - Sum of Absolute Deviations (SAD), Laplace distribution

dist1 = 3; % Error structure type:
% 0 - Normal distribution (method1 = 0)
% 1 - Poisson error structure (method1 = 0 or method1 = 1)
% 2 - Negative Binomial (VAR = factor1 * mean, empirically estimated)
% 3 - Negative Binomial (VAR = mean + alpha * mean)
% 4 - Negative Binomial (VAR = mean + alpha * mean^2)
% 5 - Negative Binomial (VAR = mean + alpha * mean^d)
% 6 - SAD method, Laplace distribution (method1 = 6)

% Automatically assign dist1 based on method1
switch method1
    case 1
        dist1 = 1;
    case 3
        dist1 = 3;
    case 4
        dist1 = 4;
    case 5
        dist1 = 5;
    case 6
        dist1 = 6;
end

numstartpoints = 10; % Number of initial guesses for optimization using MultiStart
B = 300;             % Number of bootstrap realizations for parameter uncertainty

% <============================================================================>
% <============================== ODE Model ==================================>
% <============================================================================>
% Define the model function and associated parameters for the SEIR model.

model.fc = @SEIR1;           % Name of the model function
model.name = 'SEIR model';   % Name of the ODE model

% Define model parameters:
params.label = {'\beta', '\kappa', '\gamma', 'N'}; % Symbols for model parameters
params.LB = [0.01, 0.01, 0.01, 20];              % Lower bounds for parameter estimates
params.UB = [10, 2, 2, 1000000];                 % Upper bounds for parameter estimates
params.initial = [0.6, 1/1.9, 1/4.1, 550000];    % Initial guesses for parameter values
params.fixed = [0, 1, 1, 1];                     % Boolean vector to fix parameters (1) or estimate (0)
params.fixI0 = 1;                                % Fix initial value of fitting variable to first observation (1 = yes, 0 = no)
params.composite = @R0s;                         % Function to estimate composite parameter (e.g., R_0)
params.composite_name = 'R_0';                   % Name of the composite parameter
params.extra0 = [];                              % Placeholder for additional parameters

% Define model variables:
vars.label = {'S', 'E', 'I', 'R', 'C'};             % Symbols for model variables
vars.initial = [params.initial(4) - 4, 0, 4, 0, 4]; % Initial conditions for model variables
vars.fit_index = 5;                                 % Index of the variable to fit to observed data
vars.fit_diff = 1;                                  % Fit the derivative of the fitting variable (1 = yes, 0 = no)

% <============================================================================>
% <======= Parameters for Rolling Window Analysis ===========================>
% <============================================================================>
% Settings for rolling window analysis.

windowsize1 = 17; % Size of the moving window (e.g., 17 time units)
tstart1 = 1;      % Start time point for rolling window analysis
tend1 = 1;        % End time point for rolling window analysis
printscreen1 = 1; % Boolean: Print results to screen (1 = yes, 0 = no)

end
