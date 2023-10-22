%% EXAMPLE_Amen_Galaxies_OEFPIL2D

clear
close all

%% LOAD DATA

load DATA_Amen_Galaxies;

%%  Set measurements x and y and uncertainty matrix

x = data(:,1);
y = data(:,2);

Ux  = diag(xcov);
Uy  = diag(ycov);
Uxy = diag(xycov);
U   = [Ux Uxy;Uxy Uy];

%% Set the function of parameter constraints
%  fun = fun_munubeta;

fun    =  @(mu,nu,beta) beta(1)*mu + beta(2) - nu;

%% Set the derivatives of the function of parameter constraints (optional)

options.funDiff_mu   = funDiff_mu;
options.funDiff_nu   = funDiff_nu;
options.funDiff_beta = funDiff_beta;

%% Other options

options.criterion = 'parameterdifferences';

%% Set the starting values of mu, nu and beta

mu0   = x;
nu0   = y;
beta0 = start; 

%% Run the OEFPIL algorithm

result = OEFPIL2D(x,y,U,fun,mu0,nu0,beta0,options);

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
legend('observed values','fitted values','fitted function','Location','northwest')
title('Amen Galaxies Data: Observed vs. fitted values')