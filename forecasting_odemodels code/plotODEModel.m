function [Ys,curves]=plotODEModel(options_pass,windowsize1_pass)

% <============================================================================>
% < Author: Gerardo Chowell ==================================================>
% <============================================================================>

close all
% <============================================================================>
% <=================== Declare global variables ===============================>
% <============================================================================>

global method1 % Parameter estimation method

% <============================================================================>
% <=================== Load parameter values supplied by user =================>
% <============================================================================>

if exist('options_pass','var')==1 && isempty(options_pass)==0

    options1=options_pass;

    [cadfilename1_INP,caddisease_INP,datatype_INP, dist1_INP, numstartpoints_INP,M_INP, model_INP, params_INP, vars_INP, windowsize1_INP, tstart1_INP, tend1_INP, printscreen1_INP]=options1();

else

    [cadfilename1_INP,caddisease_INP,datatype_INP, dist1_INP, numstartpoints_INP,M_INP, model_INP, params_INP, vars_INP, windowsize1_INP,tstart1_INP,tend1_INP,printscreen1_INP]=options_fit;

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
    end
end

if length(params.label)~=params.num | length(params.fixed)~=params.num | length(params.LB)~=params.num | length(params.UB)~=params.num | length(params.initial)~=params.num
    error('one or more parameter vectors do not match the number of parameters specified in <params.num>')
end


if exist('windowsize1_pass','var')==1 && isempty(windowsize1_pass)==0

    windowsize1=windowsize1_pass;
else
    windowsize1=windowsize1_INP;
end

printscreen1=printscreen1_INP;

% <==============================================================================>
% <======================== Load epidemic data ========================================>
% <==============================================================================>

[~, baseName, ext] = fileparts(cadfilename1);
if isempty(ext) || ~strcmpi(ext, '.txt')
    cadfilename1 = [baseName, '.txt'];
else
    cadfilename1 = [baseName, '.txt'];  % enforce lowercase .txt
end

%  Build platform-independent path with fullfile
fullFilePath = fullfile('.', 'input', cadfilename1);

if isfile(fullFilePath)
    % File exists.
    data=load(fullFilePath);

    data=data(1:1:windowsize1,:);

else
    warning('File does not exist. Creating a data file with zeros.')
    data=[(0:1:windowsize1-1)' zeros(windowsize1,length(vars.fit_index))];
end


figure

windowsize1

data

% solve model numerically
%timevect=tstart1:1:tstart1+windowsize1-1;
timevect=data(:,1)


options=[];
IC=vars.initial;

cc1=1;

for x=1:length(vars.fit_index)

    curvess=[];
    composite1=[];

    SSEs=[];

    tiledlayout(length(vars.fit_index),2,'Padding', 'compact', 'TileSpacing', 'compact')

    for j=1:M

        param1=[];

        for i=1:params.num

            if params.fixed(i) || j==1
                param1=[param1;params.initial(i)];

            else
                param1=[param1;unifrnd(params.LB(i),params.UB(i))];
            end

        end

        if isempty(params.composite)==1
            composite1=[composite1;NaN];
        else
            composite1=[composite1;params.composite(param1')];
        end
 
        [~,F]=ode15s(model.fc,timevect,IC,options,param1,params.extra0);

        F=real(F);

        for i2=1:vars.num
            Ys(i2,j)={F(:,i2)};
        end

        %plot(cell2mat(Ys(1,9,:))) %M,var,time

        if vars.fit_diff(x)==1
            fitcurve=abs([F(1,vars.fit_index(x));diff(F(:,vars.fit_index(x)))]);
        else
            fitcurve=F(:,vars.fit_index(x));
        end

        nexttile(cc1)
        % plot the fitting variable
        line1=plot(timevect,fitcurve,'b-');
        hold on
        nexttile(cc1+1)
        line1=plot(timevect,fitcurve,'b-');
        hold on
        curvess=[curvess fitcurve];

        % compute SSE across curves
        SSEs=[SSEs;sum((fitcurve-data(:,x+1)).^2)];

    end

    [min1,index1]=min(SSEs);

    nexttile(cc1)    % plot the fitting variable
    line1=plot(timevect,curvess(:,index1),'g-');
    set(line1,'LineWidth',4)
    hold on
    nexttile(cc1+1)
    line1=plot(timevect,curvess(:,index1),'g-');
    set(line1,'LineWidth',4)

    nexttile(cc1)
    xlabel('Time')

    if vars.fit_diff(x)
        ylabel(strcat(vars.label(vars.fit_index(x)),'''(t)',{' '}))
    else
        ylabel(strcat(vars.label(vars.fit_index(x)),'(t)',{' '}))
    end

    nexttile(cc1+1)
    xlabel('Time')
    if vars.fit_diff(x)
        ylabel(strcat(vars.label(vars.fit_index(x)),'''(t)',{' '}))
    else
        ylabel(strcat(vars.label(vars.fit_index(x)),'(t)',{' '}))
    end

    %subplot(1,2,1)
    %line1=plot(timevect,median(curvess,2),'k--')
    %set(line1,'LineWidth',3)
    %subplot(1,2,2)
    %line1=plot(timevect,median(curvess,2),'k--')
    %set(line1,'LineWidth',3)

    % plot time series data

    nexttile(cc1)
    line1=plot(data(:,1),data(:,x+1),'ro');
    set(line1,'markersize',6,'LineWidth',2)
    axis([data(1,1) data(end,1) 0 max(max(curvess))+5])

    nexttile(cc1+1)
    line1=plot(data(:,1),data(:,x+1),'ro');
    set(line1,'markersize',6,'LineWidth',2)
    axis([data(1,1) data(end,1) 0 max(data(:,x+1))+5])

    nexttile(cc1)
    set(gca,'FontSize',GetAdjustedFontSize);
    set(gcf,'color','white')

    nexttile(cc1+1)
    set(gca,'FontSize',GetAdjustedFontSize);
    set(gcf,'color','white')
    title('zoomed')

    cc1=cc1+2;

end


% plot the empirical distribution of the composite parameter

if isempty(params.composite)==0

    figure
    tiledlayout(1,1,'Padding', 'compact', 'TileSpacing', 'compact')

    composite=[median(composite1) quantile(composite1,0.025) quantile(composite1,0.975)]

    cad1=strcat('\it{',params.composite_name,'}=',num2str(composite(end,1),4),' (95%CI:',num2str(composite(end,2),4),',',num2str(composite(end,3),4),')')
    hist(composite1)
    ylabel('Frequency')
    xlabel(params.composite_name)
    title(cad1)

    set(gca,'FontSize',GetAdjustedFontSize);
    set(gcf,'color','white')

end


%% plot all state variables in a figure

figure

factor1=factor(vars.num);

if length(factor1)==1
    rows1=1;
    cols1=factor1;
elseif length(factor1)>2
    rows1=factor1(1)*factor1(2);
    cols1=factor1(3);
else
    rows1=factor1(1);
    cols1=factor1(2);
end

 tiledlayout(rows1,cols1,'Padding', 'compact', 'TileSpacing', 'compact')

cc1=1;
for i=1:1:vars.num

    nexttile(cc1)

    plot(cell2mat(Ys(i,:,:)),'b-')
    hold on

    %for j=1:M
    %    plot(cell2mat(Ys(j,i,:)),'b-')
    %    hold on
    %end

    title(vars.label(i))
    set(gca,'FontSize',GetAdjustedFontSize);
    set(gcf,'color','white')

    cc1=cc1+1;

end

for j=1:1:cols1

    nexttile(rows1*cols1-cols1+j)
    xlabel('Time')
end


if 1

    % <===============================================================================================================>
    % <=========================== Generate simulated data with a given error structure ===============================>
    % <================================================================================================================>

    %dist1=0; % Normal distribution to model error structure (method1=0)
    %dist1=1; % Poisson error structure (method1=0 OR method1=1)
    %dist1=2; % Neg. binomial error structure where var = factor1*mean where
    % factor1 is empirically estimated from the time series
    % data (method1=0)
    %dist1=3; % MLE (Neg Binomial) with VAR=mean+alpha*mean  (method1=3)
    %dist1=4; % MLE (Neg Binomial) with VAR=mean+alpha*mean^2 (method1=4)
    %dist1=5; % MLE (Neg Binomial)with VAR=mean+alpha*mean^d (method1=5)

    switch dist1

        case 0

          factor1=4; %Normal error structure

        case 3

            factor1=3;

        case 4

            factor1=3;

    end

    %dist1=6; factor1=20; %Laplace error structure

    %dist1=1;
    %factor1=1;

    d=1;

    curves=[];
    for j=1:length(vars.fit_index)
        [~,F]=ode15s(model.fc,timevect,IC,options,params.initial,params.extra0);
        if vars.fit_diff(j)==1
            fitcurve=abs([F(1,vars.fit_index(j));diff(F(:,vars.fit_index(j)))]);
        else
            fitcurve=F(:,vars.fit_index(j));
        end

        curve_noise=AddErrorStructure(cumsum(fitcurve),1,dist1,factor1,d)

        curves=[curves curve_noise];
    end

    curves=max(curves,0);

    %plot(timevect,curves)
    %hold on

    curves=[timevect curves]

    %save(strcat('./input/curve-',model.name,'-M-',num2str(M),'-dist1-',num2str(dist1),'-factor1-',num2str(factor1),'.txt'),'curves','-ascii')

end


%%%%

% Display Parameters in the Command Window
disp('<============================================================================>');
disp('                          Parameter Settings Summary                          ');
disp('<============================================================================>');

% Display Dataset Properties
disp('Dataset Properties:');
disp(['  - Time-series Data File: ', cadfilename1]);
disp(['  - Disease: ', caddisease]);
disp(['  - Data Type: ', datatype]);
disp('<============================================================================>');

% Display Parameter Estimation Settings
disp('Parameter Estimation Settings:');
disp(['  - Estimation Method (method1): ', num2str(method1)]);
switch method1
    case 0
        disp('    Method Description: Nonlinear Least Squares (LSQ)');
    case 1
        disp('    Method Description: Maximum Likelihood Estimation (MLE) Poisson');
    case 3
        disp('    Method Description: MLE Negative Binomial (VAR = mean + alpha * mean)');
    case 4
        disp('    Method Description: MLE Negative Binomial (VAR = mean + alpha * mean^2)');
    case 5
        disp('    Method Description: MLE Negative Binomial (VAR = mean + alpha * mean^d)');
    case 6
        disp('    Method Description: Sum of Absolute Deviations (SAD), Laplace distribution');
    otherwise
        disp('    Method Description: Unknown');
end

disp(['  - Error Structure (dist1): ', num2str(dist1)]);
switch dist1
    case 0
        disp('    Error Structure Description: Normal Distribution');
    case 1
        disp('    Error Structure Description: Poisson Error Structure');
    case 2
        disp('    Error Structure Description: Negative Binomial (VAR = factor1 * mean)');
    case 3
        disp('    Error Structure Description: Negative Binomial (VAR = mean + alpha * mean)');
    case 4
        disp('    Error Structure Description: Negative Binomial (VAR = mean + alpha * mean^2)');
    case 5
        disp('    Error Structure Description: Negative Binomial (VAR = mean + alpha * mean^d)');
    case 6
        disp('    Error Structure Description: Laplace Distribution (SAD)');
    otherwise
        disp('    Error Structure Description: Unknown');
end

disp(['  - Number of Initial Guesses (MultiStart): ', num2str(numstartpoints)]);
disp(['  - Number of Bootstrap Realizations: ', num2str(M)]);
disp('<============================================================================>');

% Display ODE Model Information
disp('ODE Model:');
disp(['  - Model Function: ', func2str(model.fc)]);
disp(['  - Model Name: ', model.name]);
disp('<============================================================================>');

% Display Model Parameters
disp('Model Parameters:');
disp(['  - Labels: ', strjoin(params.label, ', ')]);
disp(['  - Lower Bounds: ', num2str(params.LB)]);
disp(['  - Upper Bounds: ', num2str(params.UB)]);
disp(['  - Initial Guesses: ', num2str(params.initial)]);
disp(['  - Fixed Parameters: ', num2str(params.fixed)]);
disp(['  - Fix Initial Value (fixI0): ', logicalToString(params.fixI0)]);
disp('<============================================================================>');

% Display Model Variables
disp('Model Variables:');
disp(['  - Labels: ', strjoin(vars.label, ', ')]);
disp(['  - Initial Conditions: ', num2str(vars.initial)]);
disp(['  - Fit Variable Index: ', num2str(vars.fit_index)]);
disp(['  - Fit Derivative (fit_diff): ', logicalToString(vars.fit_diff)]);
disp('<============================================================================>');

% Display Rolling Window Parameters
disp('Rolling Window Parameters:');
disp(['  - Window Size: ', num2str(windowsize1)]);
disp(['  - Start Time Point: ', num2str(tstart1_INP)]);
disp(['  - End Time Point: ', num2str(tend1_INP)]);
disp(['  - Print Results to Screen (printscreen1): ', logicalToString(printscreen1)]);
disp('<============================================================================>');
disp('                          End of Parameter Summary                            ');
disp('<============================================================================>');

% Logical to String Conversion
    function str = logicalToString(logicalValue)
        if logicalValue
            str = 'Yes';
        else
            str = 'No';
        end
    end

end

