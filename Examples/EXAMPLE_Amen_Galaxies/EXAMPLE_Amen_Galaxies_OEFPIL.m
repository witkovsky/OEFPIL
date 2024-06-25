%% EXAMPLE_Amen_Galaxies_OEFPIL

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
%  fun = fun_mubeta;

fun    =  @(mu,beta) beta(1)*mu{1} + beta(2) - mu{2};

%% Set the derivatives of the function of parameter constraints (optional)

options.funDiff_mu   = @(mu,beta) {beta(1)*ones(size(mu{1})), -ones(size(mu{2}))};
options.funDiff_beta = @(mu,beta) [mu{1},ones(size(mu{1}))];

%% Other options

options.criterion = 'parameterdifferences';

%% Set the starting values of mu, nu and beta

mu0   = {x, y};
beta0 = start; 

%% Run the OEFPIL algorithm

result = OEFPIL({x,y},U,fun,mu0,beta0,options);

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
title('Amen Galaxies Data: Observed vs. fitted values')