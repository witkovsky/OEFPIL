%% EXAMPLE_FlowMeter / Fractional Polynomial Calibration

clear
close all

%% LOAD DATA

load DATA_FlowMeter;

%%  Set measurements x and y and uncertainty matrix

x = data(:,1);
y = data(:,2);

Ux  = xcov;
Uy  = ycov;
Uxy = [];

%% Set the function of parameter constraints
%  fun = fun_munubeta;

fun  = @(mu,nu,beta) beta(1)*mu.^(-1) + beta(2)*mu.^(-0.5) ...
       + beta(3) + beta(4)*mu.^(0.5) + beta(5)*mu - nu;

%% Set the derivatives of the function of parameter constraints (optional)

options.funDiff_mu   = funDiff_mu;
options.funDiff_nu   = funDiff_nu;
options.funDiff_beta = funDiff_beta;

%% Other options

options.tol = 1e-10;
options.method = 'oefpil';
options.criterion = 'parameterdifferences';

%% Set the starting values of mu, nu and beta

mu0   = x;
nu0   = y;
beta0 = start; 

%% Run the OEFPIL algorithm

result = OEFPIL2D(x,y,{Ux,Uy,Uxy},fun,mu0,nu0,beta0,options);

%% Plot Observed data + Fitted values + Fitted line

xx = linspace(min(x),max(x),51)';
yy = fun(xx,xx,result.beta) + xx;

figure
plot(x,y,'*')
hold on
plot(result.mu,result.nu,'o')
plot(xx,yy,'-')
grid on
xlabel('x')
ylabel('y')
legend('observed values','fitted values','fitted function','Location','northeast')
title('Flow Meter Data: Observed vs. fitted values')
