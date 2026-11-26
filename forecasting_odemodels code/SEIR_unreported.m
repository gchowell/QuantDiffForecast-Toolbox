% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

function dx=SEIR_unreported(t,x,params0,extra0)

beta0=params0(1);
beta1=params0(2);
alpha=params0(3); % non-homogenous mixing parameter
q1=params0(4);
rho=params0(5);
k=params0(6);
gamma1=params0(7);
N=params0(8);

t0=35;

if t<t0

    betaf=beta0;

else

    betaf=beta1+(beta0-beta1)*exp(-q1*(t-t0)); % exponential decline

    %betaf=beta0*exp(-q1*t); % exponential decline
    %betaf=beta0*(1./(t+1)).^q1; % power-law decline

end

dx=zeros(6,1);  % define the vector of the state derivatives: S, E, I, U, R, C

dx(1,1)= -betaf*x(1,1).*((x(3,1)+x(4,1)).^alpha)./N; %S

dx(2,1)= betaf*x(1,1).*((x(3,1)+x(4,1)).^alpha)./N - k*x(2,1); %E

dx(3,1)= k*rho*x(2,1) - gamma1*x(3,1); %I

dx(4,1)= k*(1-rho)*x(2,1) - gamma1*x(4,1); %U

dx(5,1)= gamma1*(x(3,1)+x(4,1)); %R

dx(6,1)= k*rho*x(2,1); %C

