%% EXAMPLE_OliverPharr_OEFPIL2D / Function fit fun(mu,nu,[0.75;-0.25;1.75])
%  Artificial generated data

clear
close all

%% LOAD DATA

load DATA_OliverPharr;

%%  Set measurements x and y and uncertainty matrix

 x     = [ 0.2505    2.6846    5.1221    7.5628   10.0018 ...
           0.2565    2.6858    5.1255    7.5623    9.9952 ...
           0.2489    2.6830    5.1271    7.5603   10.0003]';

 y     = [ 0.2398    4.9412   14.2090   27.3720   44.0513 ...
           0.2110    4.9517   14.2306   27.3937   44.0172 ...
           0.2303    4.9406   14.2690   27.3982   44.0611]';

%% Set the function of parameter constraints

 fun   = @(mu,nu,beta) beta(1).*(mu-beta(2)).^beta(3) - nu;
 
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
legend('observed values','fitted values','fitted function','Location','northwest')
title('Oliver Pharr Data: Observed vs. fitted values')
