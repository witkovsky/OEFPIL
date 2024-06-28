%% EXAMPLE_FlowMeter / Fractional Polynomial Calibration

clear
close all

%% LOAD DATA

load DATA_FlowMeter;

%%  Set measurements x and y and uncertainty matrix

x = data(:,1);
y = data(:,2);
data = {x, y};

Ux  = xcov;
Uy  = ycov;
U   = {Ux []; [] Uy};

%% Set the function of parameter constraints

fun  = @(mu,beta) beta(1)*mu{1}.^(-1) + beta(2)*mu{1}.^(-0.5) ...
       + beta(3) + beta(4)*mu{1}.^(0.5) + beta(5)*mu{1} - mu{2};

%% Set the derivatives of the function of parameter constraints (optional)

options.funDiff_mu   = @(mu,beta) {-beta(1)*mu{1}.^(-2) ...
    - 0.5*beta(2)*mu{1}.^(-1.5) + 0.5*beta(4)*mu{1}.^(-0.5) ...
    + beta(5)*ones(size(mu{1})), -ones(size(mu{2}))};
options.funDiff_beta = @(mu,beta) [mu{1}.^(-1),mu{1}.^(-0.5), ...
    ones(size(mu{1})),mu{1}.^(0.5),mu{1}];

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
legend('observed values','fitted values','fitted function','Location','northeast')
title('Flow Meter Data: Observed vs. fitted values')