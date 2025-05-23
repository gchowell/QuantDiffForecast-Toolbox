
function Run_PracticalIndentifiability_ODEModel(options_pass,factor1,replicates1)

%Run_PracticalIndentifiability_ODEModel(@options_forecast_SEIR_plain_unreported_dist1_1,4,2)

% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

% Fitting and forecasting model to epidemic data with quantified uncertainty

close all

% <============================================================================>
% <=================== Declare global variables ===============================>
% <============================================================================>

global method1 % Parameter estimation method

% <============================================================================>
% <=================== Load parameter values supplied by user =================>
% <============================================================================>

if exist('options_pass','var')==1 && isempty(options_pass)==0

    options=options_pass; %forecast horizon (number of data points ahead)

    [cadfilename1_INP,caddisease_INP,datatype_INP, dist1_INP, numstartpoints_INP,M_INP, model_INP, params_INP, vars_INP, getperformance_INP,forecastingperiod_INP,windowsize1_INP,tstart1_INP,tend1_INP,printscreen1_INP]=options();

else

    options=str2func('options_forecast.m');

    [cadfilename1_INP,caddisease_INP,datatype_INP, dist1_INP, numstartpoints_INP,M_INP, model_INP, params_INP, vars_INP, getperformance_INP,forecastingperiod_INP,windowsize1_INP,tstart1_INP,tend1_INP,printscreen1_INP]=options_forecast;

end

params_INP.num=length(params_INP.label); % number of model parameters

vars_INP.num=length(vars_INP.label); % number of variables comprising the ODE model


% <============================================================================>
% <================================ Datasets properties ==============================>
% <============================================================================>

cadfilename1=cadfilename1_INP;

DT=1;

caddisease=caddisease_INP;

datatype=datatype_INP;

% <=============================================================================>
% <=========================== Parameter estimation ============================>
% <=============================================================================>

%method1=0; % Type of estimation method: 0 = LSQ

d=1;

dist1=dist1_INP; %Define dist1 which is the type of error structure:

if method1>0
    dist1=method1;
end

% LSQ=0,
% MLE Poisson=1,
% Pearson chi-squard=2,
% MLE (Neg Binomial)=3, with VAR=mean+alpha*mean;
% MLE (Neg Binomial)=4, with VAR=mean+alpha*mean^2;
% MLE (Neg Binomial)=5, with VAR=mean+alpha*mean^d;


numstartpoints=numstartpoints_INP; % Number of initial guesses for optimization procedure using MultiStart

M=M_INP; % number of bootstrap realizations to characterize parameter uncertainty


% <==============================================================================>
% <============================== ODE model =====================================>
% <==============================================================================>

model=model_INP;
params=params_INP;
vars=vars_INP;

for j=1:params.num
    if params.initial(j)<params.LB(j) | params.initial(j)>params.UB(j)
        error('values in <params.initial> should lie within their parameter bounds defined by <params.LB> and <params.UB> ')
        return
    end
end

if length(params.label)~=params.num | length(params.fixed)~=params.num | length(params.LB)~=params.num | length(params.UB)~=params.num | length(params.initial)~=params.num
    error('one or more parameter vectors do not match the number of parameters specified in <params.num>')
end


% <==============================================================================>
% <========================== Forecasting parameters ===================================>
% <==============================================================================>

getperformance=getperformance_INP; % flag or indicator variable (1/0) to calculate forecasting performance or not

if exist('forecastingperiod_pass','var')==1 & isempty(forecastingperiod_pass)==0

    forecastingperiod=forecastingperiod_pass; %forecast horizon (number of data points ahead)

else

    forecastingperiod=forecastingperiod_INP;

end

% <==================================================================================>
% <========================== Parameters of the rolling window analysis =========================>
% <==================================================================================>

tstart1=tstart1_INP;
tend1=tend1_INP;
windowsize1=windowsize1_INP;
printscreen1=printscreen1_INP;


for j=1:replicates1

    [Ys,curves]=Run_simulate_ODEModel(options,factor1);

    vars.fit_index

    simulation=[];

    for i=1:length(vars.fit_index)

        data=cell2mat(Ys(vars.fit_index(i),1,:));

        if vars.fit_diff(i)==1
            simulation=[simulation abs([data(1);diff(data)])];
        else
            simulation=[simulation data];
        end

    end

    simulation=[(0:1:length(simulation(:,1))-1)' simulation];

    figure(200)
    plot(curves(:,1),curves(:,2:end),'bo')
    hold on
    plot(simulation(:,1),simulation(:,2:end),'r-')


    % Check if fileName ends with '.txt'
    if ~endsWith(cadfilename1, '.txt', 'IgnoreCase', true)
        % Append '.txt' extension if not present
        cadfilename1 = strcat(cadfilename1, '.txt');
    end

    % Create full file path
    fullFilePath = fullfile('./input', cadfilename1);

    save(fullFilePath,'curves','-ascii')

    % Generate forecast
    Run_Forecasting_ODEModel(options,1,1,windowsize1,forecastingperiod);

    %Read forecast and quantify change in forecast uncertainty
    for j=1:length(vars.fit_index)

        T=readtable(strcat('./output/Forecast-model_name-',model.name,'-vars.fit_index-',num2str(vars.fit_index(j)),'-tstart-',num2str(tstart1),'-fixI0-',num2str(params.fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-horizon-',num2str(forecastingperiod),'-',caddisease,'-',datatype,'.csv'));

        width_calib=T.UB(1:windowsize1)-T.LB(1:windowsize1)

        width_pred=T.UB(windowsize1+1:end)-T.LB(windowsize1+1:end)

        % Find estimated horizon at which the width of the 95%PI doubles
        % starting from the last data point of the calibration period
        estimated_horizon = find_time_for_value(1:forecastingperiod,width_pred, 2*width_calib(end));

        % Display result
        fprintf('Estimated horizon for value %.2f: %.2f\n', 2*width_calib(end), estimated_horizon);

    end

end

