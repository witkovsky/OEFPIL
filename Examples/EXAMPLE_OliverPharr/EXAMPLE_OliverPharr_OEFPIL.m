%% EXAMPLE_OliverPharr_OEFPIL / Function fit fun(mu,nu,[0.75;-0.25;1.75])
%  Artificial generated data

clear
close all

%% LOAD DATA

load DATA_OliverPharr.mat;

%%  Set measurements x and y and uncertainty matrix

x     = [ 0.2505    2.6846    5.1221    7.5628   10.0018 ...
          0.2565    2.6858    5.1255    7.5623    9.9952 ...
          0.2489    2.6830    5.1271    7.5603   10.0003]';

y     = [ 0.2398    4.9412   14.2090   27.3720   44.0513 ...
          0.2110    4.9517   14.2306   27.3937   44.0172 ...
          0.2303    4.9406   14.2690   27.3982   44.0611]';

data = {x, y};

U   = {Ux []; [] Uy};

%% Set the function of parameter constraints

fun = @(mu,beta) beta(1).*(mu{1}-beta(2)).^beta(3) - mu{2};
 
%% Set the derivatives of the function of parameter constraints (optional)

options.funDiff_mu = @(mu,beta) { ...
    beta(1)*beta(3)*(mu{1}-beta(2)).^(beta(3)-1), ...
    -ones(size(mu{1}))};
options.funDiff_beta = @(mu,beta) [(mu{1}-beta(2)).^beta(3), ...
    -beta(1)*beta(3)*(mu{1}-beta(2)).^(beta(3)-1), ...
    beta(1).*(mu{1}-beta(2)).^beta(3).*log(mu{1}-beta(2))];

%% Other options

options.tol = 1e-10;
options.method = 'oefpil';
options.criterion = 'parameterdifferences';

%% Set the starting values of mu, nu and beta

mu0   = {x, y};
beta0 = start; 
 
%% Run the OEFPIL algorithm

result = OEFPIL(data,U,fun,mu0,beta0,options);

%% Plot Observed data + Fitted values + Fitted line

xx = linspace(min(x),max(x),51)';
yy = fun({xx,xx},result.beta) + xx;

figure
plot(x,y,'*')
hold on
plot(result.muCell{1},result.muCell{2},'o')
plot(xx,yy,'-')
grid on
xlabel('x')
ylabel('y')
legend('observed values','fitted values','fitted function','Location','northwest')
title('Oliver Pharr Data: Observed vs. fitted values')