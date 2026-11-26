% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

function [cadfilename1,caddisease,datatype, dist1, numstartpoints,B, model, params,vars,windowsize1,tstart1,tend1, printscreen1]=options_fit_SEIR1


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

% This example fits SIMULATED data, so the true value of the estimated parameter
% is known: 100 daily counts from SEIR1 with beta=0.6, kappa=1/2, gamma=1/3,
% N=100000 and initial state [N-5 0 5 0 5] (the values in params.initial and
% vars.initial below), with Poisson noise on the daily increments of C after the
% first point. A fit should return beta close to 0.6 (R_0 = beta/gamma = 1.8).
% The series peaks near row 66; the default window uses the first 17 rows, where
% counts are small, so try a longer window as well, e.g.
%   Run_Fit_ODEModel(@options_fit_SEIR1, 1, 1, 40)

cadfilename1='curve-SEIR-simulated-poisson'; % Name of the data file containing the incidence curve

caddisease='simulated'; % string indicating the name of the disease related to the time series data

datatype='cases'; % string indicating the nature of the data (cases, deaths, hospitalizations, etc)

% <=============================================================================>
% <=========================== Parameter estimation ============================>
% <=============================================================================>

method1=1; % Type of estimation method

% Nonlinear least squares (LSQ)=0,
% MLE Poisson=1,
% MLE (Neg Binomial)=3, with VAR=mean+alpha*mean;
% MLE (Neg Binomial)=4, with VAR=mean+alpha*mean^2;
% MLE (Neg Binomial)=5, with VAR=mean+alpha*mean^d;

dist1=1; % Define dist1 which is the type of error structure. See below:

%dist1=0; % Normal distribution to model error structure (method1=0)
%dist1=1; % Poisson error structure (method1=0 OR method1=1)
%dist1=2; % Neg. binomial error structure where var = factor1*mean where
                  % factor1 is empirically estimated from the time series
                  % data (method1=0)
%dist1=3; % MLE (Neg Binomial) with VAR=mean+alpha*mean  (method1=3)
%dist1=4; % MLE (Neg Binomial) with VAR=mean+alpha*mean^2 (method1=4)
%dist1=5; % MLE (Neg Binomial)with VAR=mean+alpha*mean^d (method1=5)


switch method1
    case 1
        dist1=1;
    case 3
        dist1=3;
    case 4
        dist1=4;
    case 5
        dist1=5;
end


numstartpoints=10; % Number of initial guesses for optimization procedure using MultiStart

B=100; % number of bootstrap realizations to characterize parameter uncertainty

% <==============================================================================>
% <============================== ODE model =====================================>
% <==============================================================================>

model.fc=@SEIR1; % name of the model function
model.name='SEIR model';   % string indicating the name of the ODE model

params.num=4; % number of model parameters
params.label={'\beta','\kappa','\gamma','N'};  % list of symbols to refer to the model parameters
params.LB=[0.01 0.01 0.01 20]; % lower bound values of the parameter estimates
params.UB=[10 2 2 1000000]; % upper bound values of the parameter estimates
params.initial=[0.6 1/2 1/3 100000]; % initial parameter values/guesses
params.fixed=[0 1 1 1]; % Boolean vector to indicate any parameters that should remain fixed (1) to initial values indicated in params.initial. Otherwise the parameter is estimated (0).
params.fixI0=1; % Boolean variable indicating if the initial value of the fitting variable is fixed according to the first observation in the time series (1). Otherwise, it will be estimated along with other parameters.
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

vars.num=5; % number of variables comprising the ODE model
vars.label={'S','E','I','R','C'}; % list of symbols to refer to the variables included in the model
vars.initial=[params.initial(4)-5 0 5 0 5];  % vector of initial conditions for the model variables
vars.fit_index=5; % index of the model's variable that will be fit to the observed time series data
vars.fit_diff=1; % boolean variable to indicate if the derivative of model's fitting variable should be fit to data.

% <==================================================================================>
% <========================== Parameters of the rolling window analysis =========================>
% <==================================================================================>

windowsize1=17;  % moving window size

tstart1=1; % time point for the start of rolling window analysis

tend1=1;  %time point for the end of the rolling window analysis

printscreen1=0;
