% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

function [P, residual, fitcurve, forecastcurve, timevect2,initialguess,fval,F1,F2]=fit_model(data1,params0,numstartpoints,DT,modelX,paramsX,varsX,forecastingperiod)

global model params vars method1 timevect ydata

model=modelX;
params=paramsX;
vars=varsX;

% Calculate time vector from data column and a scaling factor DT
timevect = data1(:,1) * DT;

% Extract initial conditions from the second column onwards of the first row
I0 = data1(1, 2:end);

% Copy initial parameter values
z = params0;

% Adjust lower and upper bounds for parameters that are fixed
for i = 1:params.num
    if params.fixed(i)
        params.LB(i) = params0(i);
        params.UB(i) = params0(i);
    end
end

% Set extended bounds based on the method specified
switch method1
    case {0, 1,6}
        % Cases 0 and 1 have no extension on bounds
        LBe = [0 0];
        UBe = [0 0];
    case {3, 4}
        % Cases 3 and 4 allow wide variation
        LBe = [1e-8, 1];
        UBe = [1e4, 1];
    case 5
        % Case 5 has specific limits for d, assuming it should be >=1
        LBe = [1e-8, 0.6];
        UBe = [1e4, 1e3];
end

% Configure bounds based on whether initial conditions are fixed
if params.fixI0 == 1
    LB = [params.LB, I0, LBe];
    UB = [params.UB, I0, UBe];
else
    % If not fixed, use zero and sum of absolute values for initial conditions
    LB = [params.LB, zeros(1, length(I0))+0.001, LBe];
    UB = [params.UB, sum(abs(data1(:, 2:end))), UBe];
end


% 
% if 0 % USE LSQCURVEFIT (Non-linear least squares)
% 
%     options=optimset('tolfun',10^-5,'TolX',10^-5,'MaxFunEvals',3200,'MaxIter',3200, 'algorithm','trust-region-reflective');
% 
%     [P,resnorm,residual,exitflag,output,lambda,J]=lsqcurvefit(@plotModifiedLogisticGrowth1,z,timevect,data1(:,2),LB,UB,options);
% 
%     f=@plotModifiedLogisticGrowth1;
% 
%     problem = createOptimProblem('lsqcurvefit','x0',z,'objective',f,'lb',LB,'ub',UB,'xdata',timevect,'ydata',data1(:,2),'options',options);
% 
%     %ms = MultiStart('PlotFcns',@gsplotbestf,'Display','final');
% 
%     ms = MultiStart('Display','final');
% 
%     ms = MultiStart(ms,'StartPointsToRun','bounds')
% 
%     [P,errormulti] = run(ms,problem,20)
% 
%     z=P;
% 
%     [P,resnorm,residual,exitflag,output,lambda,J]=lsqcurvefit(@plotModifiedLogisticGrowth1,z,timevect,data1(:,2),LB,UB,options);
% 
% 
% end

%A=[];      % We are using fmincon, but using none of the constraint options
%b=[];
%Aeq=[];
%beq=[];
%nonlcon=[];

%options=optimset('tolfun',10^-5,'TolX',10^-5,'MaxFunEvals',3200,'MaxIter',3200, 'algorithm','interior-point');

%[P, fval, exitflag]=fmincon(@plotModifiedLogisticGrowthMethods1,z,A,b,Aeq,beq,LB,UB,nonlcon,options);

%method1=1; %LSQ=0, MLE (Poisson)=1, Pearson chi-squared=2. MLE(neg binomial)=3

ydata=data1(:,2:end);

%convert a matrix into a long column vector using the (:) 
if length(I0)>1
    ydata=ydata(:);
end

% Define the objective function handle
f = @parameterSearchODE;

% Set optimization options for fmincon
options = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...  % sqp, interior-point
    'StepTolerance', 1e-4, ... % Slightly relaxed
    'FunctionTolerance', 1e-4, ... 
    'OptimalityTolerance', 1e-4, ...
    'MaxFunEvals', 10000, ...
    'MaxIter', 10000);

% Define the optimization problem
problem = createOptimProblem('fmincon', 'objective', f, 'x0', z, ...
                             'lb', LB, 'ub', UB, 'options', options);

% Setup MultiStart with no output display
ms = MultiStart('Display', 'off');

% Define start point sets
tpoints = CustomStartPointSet(z);

% Initialize the flag to manage the while loop
flagg = -1;

% Run optimization until a stopping condition is met
while flagg < 0
    initialguess = [];

    % Add random start points if specified
    if numstartpoints > 0
        rpoints = RandomStartPointSet('NumStartPoints', numstartpoints);
        allpts = {rpoints, tpoints};
        initialguess = [list(rpoints, problem); z];
    else
        allpts = {tpoints};
        initialguess = [initialguess; z];
    end

    % Run the MultiStart solver
    [P, fval, flagg, outpt, allmins] = run(ms, problem, allpts);

end

% P is the vector with the estimated parameters

% Initialize the options for the ODE solver (if any specific options needed, define here)
options = [];

% Set initial conditions based on params.fixI0 flag
IC = vars.initial;
if params.fixI0 == 1
    IC(vars.fit_index) = I0;  % Fix the initial conditions to I0 for specified indices
else
    % If not fixed, use parameter values following the first 'num' parameters
    IC(vars.fit_index) = P(params.num + 1 : params.num + length(I0));
end

% Solve the differential equations using ode15s
[~, F] = ode15s(model.fc, timevect, IC, options, P, params.extra0);
F1 = F;  % Store the result in F1 (seems redundant unless F1 is used differently not shown here)

% Initialize variables for fitting the data
yfit = zeros(length(ydata), 1);
currentEnd = 0;

% Loop through each variable index to be fitted
for j = 1:length(vars.fit_index)
    % Check if differentiation is needed
    if vars.fit_diff(j) == 1
        % Compute absolute value of the derivative of the model output
        fitcurve = abs([F(1, vars.fit_index(j)); diff(F(:, vars.fit_index(j)))]);
    else
        % Use the model output directly
        fitcurve = F(:, vars.fit_index(j));
    end

    % Aggregate the fit results
    yfit(currentEnd + 1 : currentEnd + length(fitcurve)) = fitcurve;
    currentEnd = currentEnd + length(fitcurve);
end

% Final aggregated fit curve
fitcurve = yfit;

% 
% figure
% subplot(1,2,1)
% plot(ydata,'ko')
% hold on
% plot(fitcurve,'r')


% Calculate residuals
residual = fitcurve - ydata;

% Decide action based on forecasting period
if forecastingperiod < 1
    % If no extended forecasting is needed, adjust curve based on residuals
    forecastcurve = residual + ydata;
    timevect2 = timevect;  % Use the original time vector
    F2 = F1;               % Reuse the previously calculated ODE results
else
    % Define the extended time vector for forecasting
    timevect2 = data1(1,1) : DT : (data1(end,1) + forecastingperiod * DT);

    % Solve the ODE to get new results over the extended time period
    [~, F2] = ode15s(model.fc, timevect2, IC, options, P, params.extra0);

    % Initialize forecast storage based on the number of fit indices and time points
    yforecast = zeros(length(vars.fit_index) * length(timevect2), 1);
    currentEnd = 0;

    % Process each variable index specified for fitting
    for j = 1:length(vars.fit_index)
        % Determine if the variable requires differencing
        if vars.fit_diff(j) == 1
            % Calculate the absolute derivative of the forecast
            forecastcurve = abs([F2(1, vars.fit_index(j)); diff(F2(:, vars.fit_index(j)))]);
        else
            % Use the ODE output directly for the forecast
            forecastcurve = F2(:, vars.fit_index(j));
        end

        % Store the processed forecast in the yforecast array
        yforecast(currentEnd + 1 : currentEnd + length(forecastcurve)) = forecastcurve;
        currentEnd = currentEnd + length(forecastcurve);
    end

    % Assign the compiled forecast data as the final forecast curve
    forecastcurve = yforecast;
end

