function objfunction=parameterSearchODE(z)

global model params vars method1 timevect ydata

numFitIndices = length(vars.fit_index);  % Cache the number of fitting indices for clarity and efficiency

I0=z(params.num+1:params.num+numFitIndices);
alpha=z(params.num+numFitIndices+1);
d=z(params.num+numFitIndices+2);


% Initialize initial conditions from predefined variables
IC = vars.initial;

% Set specific initial conditions for indices specified in vars.fit_index
IC(vars.fit_index) = I0;

% Tight ODE tolerances keep the integration error far below the
% optimizer's stopping tolerances, so the objective stays smooth enough
% for finite-difference gradients. 'NonNegative' keeps compartments from
% drifting negative at extreme parameter draws, which otherwise produces
% complex logarithms downstream.
% NOTE: remove the 'NonNegative' line if your model has
% negative-valued state variables.
opts_ode = odeset('RelTol',1e-6, ...
                  'AbsTol',1e-8, ...
                  'NonNegative',1:length(IC));

% Guard the integration: ode15s returns a truncated solution (fewer rows
% than timevect) when it fails mid-interval, and can error outright for
% pathological parameter sets. Both cases are mapped to a large-but-finite
% penalty so the optimizer simply moves away from the bad region.
try
    [~, x] = ode15s(model.fc, timevect, IC, opts_ode, z, params.extra0);
catch
    objfunction = 1e10;
    return
end

if size(x,1) ~= length(timevect) || any(~isfinite(x(:))) || ~isreal(x)
    objfunction = 1e10;
    return
end

% Initialize yfit array to store fit results
yfit = zeros(length(ydata), 1);

% Initialize an index to keep track of the end position in yfit
currentEnd = 0;

% Loop through each fit index to process the ODE solution
for j = 1:length(vars.fit_index)
    % Check if differentiation is needed based on vars.fit_diff
    if vars.fit_diff(j) == 1
        % Calculate the absolute derivative of the ODE solution for this variable
        fitcurve = abs([x(1, vars.fit_index(j)); diff(x(:, vars.fit_index(j)))]);
    else
        % Use the ODE solution directly
        fitcurve = x(:, vars.fit_index(j));
    end

    % Append the fitcurve to the yfit array, updating currentEnd to the new position
    yfit(currentEnd + 1 : currentEnd + length(fitcurve)) = fitcurve;
    currentEnd = currentEnd + length(fitcurve);
end


eps=0.001;

%%MLE expression
%This is the negative log likelihood, name is legacy from least squares code
%Note that a term that is not a function of the params has been excluded so to get the actual
%negative log-likliehood value you would add: sum(log(factorial(sum(casedata,2))))

if sum(yfit)==0
    objfunction=10^10;%inf;
else
    %    z
    yfit(yfit==0)=eps; %set zeros to eps to allow calculation below.  Shouldn't affect solution, just keep algorithm going.


% Defend the likelihood computation against invalid data: MATLAB's
% gammaln throws an error for negative or otherwise invalid input rather
% than returning NaN, and bootstrap replicate data can contain NaN or
% negative values (e.g., normal/Laplace error structures draw negatives).
% An exception here would kill the local solver run silently inside
% MultiStart. Three layers of defense: (a) return the penalty if ydata is
% not finite; (b) use yg = max(ydata,0) inside the gammaln terms, so
% observations y <= 0 contribute zero to the gamma-ratio term
% (gammaln(0+m)-gammaln(m) = 0); (c) wrap the whole likelihood switch in
% try/catch so no data pathology can throw out of the objective.
if any(~isfinite(ydata))
    objfunction = 1e10;
    return
end
yg = max(ydata, 0);   % sanitized counts for gammaln arguments only

try
    switch method1

        case 0  %Least squares

            objfunction=sum((ydata-yfit).^2);


        case 1 % MLE for Poisson distribution (negative log-likelihood)

            objfunction=-sum(ydata.*log(yfit)-yfit);


        case 3  % MLE Negative binomial (negative log-likelihood) where sigma^2=mean+alpha*mean;

            % Vectorized log-likelihood using the identity
            %   sum_{j=0}^{y-1} log(j+m) = gammaln(y+m) - gammaln(m),
            % with m = yfit/alpha.
            m = yfit./alpha;
            sum1 = sum( gammaln(yg+m) - gammaln(m) ...
                        + ydata.*log(alpha) ...
                        - (ydata + m).*log(1+alpha) );

            objfunction=-sum1;

        case 4
            % MLE Negative binomial (negative log-likelihood) where sigma^2=mean+alpha*mean^2;

            % Vectorized via the same gammaln identity, with constant m = 1/alpha.
            m = 1/alpha;
            sum1 = sum( gammaln(yg+m) - gammaln(m) ...
                        + ydata.*log(alpha.*yfit) ...
                        - (ydata + m).*log(1+alpha.*yfit) );

            objfunction=-sum1;

        case 5
            % MLE Negative binomial (negative log-likelihood) where sigma^2=mean+alpha*mean^d;

            % Vectorized via the same gammaln identity, with
            % m = (1/alpha)*yfit.^(2-d).
            m = (1./alpha).*yfit.^(2-d);
            sum1 = sum( gammaln(yg+m) - gammaln(m) ...
                        + ydata.*log(alpha.*yfit.^(d-1)) ...
                        - (ydata + m).*log(1+alpha.*yfit.^(d-1)) );

            objfunction=-sum1;

        case 6

           objfunction=sum(abs(ydata-yfit));


    end

    % Any unexpected error in the likelihood computation is converted into
    % the standard penalty instead of terminating the local solver run.
    catch
        objfunction = 1e10;
        return
    end

end

% Never hand NaN/Inf/complex back to fmincon: line search and finite
% differencing degrade badly on non-finite values, so a large finite
% penalty makes such points cleanly rejectable.
if ~isreal(objfunction) || ~isfinite(objfunction)
    objfunction = 1e10;
end
