%% EXAMPLE_Area_Calibration_BV52_OEFPIL

clear
close all

%% LOAD DATA

load DATA_Area_Calibration_BV52

%%  Set measurements x and y and uncertainty matrix

x = hcApunc(:,1);
y = hcApunc(:,2);

U = {Ux Uxy; Uxy' Uy};


%% Set fun1 : Fractional polynomial "a2*x^2 + a1*x + a12*sqrt(x)"
% fun1 =  @(mu,beta) beta(1)*mu{1}.^2 + beta(2)*mu{1} ...
%     + beta(3)*mu{1}.^(1/2) - mu{2};
% 
% fun1Diff_mu = @(mu,beta) {2*beta(1)*mu{1} + beta(2) ...
%     + (1/2)*beta(3)*mu{1}.^(-1/2), -ones(size(mu{2}))};
% 
% fun1Diff_beta = @(mu,beta) [mu{1}.^2, mu{1}, mu{1}.^(1/2)];

%% Set fun2 : Fractional polynomial "a2*x^2 + a1*x + a12*x^(1/2) + a14*x^(1/4)"

% fun2 =  @(mu,beta) beta(1)*mu{1}.^2 + beta(2)*mu{1} ...
%     + beta(3)*mu{1}.^(1/2) + beta(4)*mu{1}.^(1/4)  - mu{2};
% 
% fun2Diff_mu = @(mu,beta) {2*beta(1)*mu{1} + beta(2) ...
%     + (1/2)*beta(3)*mu{1}.^(-1/2) + (1/4)*beta(4)*mu{1}.^(-3/4), ...
%     -ones(size(mu{2}))};
% 
% fun2Diff_beta = @(mu,beta) [mu{1}.^2, mu{1}, mu{1}.^(1/2), mu{1}.^(1/4)];

%% Set fun3 : Polynomial "a3*x^3 + a2*x^2 + a1*x"

% fun3 =  @(mu,beta) beta(1)*mu{1}.^3 + beta(2)*mu{1}.^2 + beta(3)*mu{1} ...
%         - mu{2};
% 
% fun3Diff_mu = @(mu,beta) {3*beta(1)*mu{1}.^2 + 2*beta(2)*mu{1} ...
%     + beta(3), -ones(size(mu{2}))};
% 
% fun3Diff_beta = @(mu,beta) [mu{1}.^3, mu{1}.^2, mu{1}];


%% Set fun4 : Polynomial "a4*x^4 + a3*x^3 + a2*x^2 + a1*x"

% fun4 =  @(mu,beta) beta(1)*mu{1}.^4 + beta(2)*mu{1}.^3 + ...
%          + beta(3)*mu{1}.^2 + beta(4)*mu{1} - mu{2};
% 
% fun4Diff_mu = @(mu,beta) {4*beta(1)*mu{1}.^3 + 3*beta(2)*mu{1}.^2 ...
%     + 2*beta(3)*mu{1} + beta(4), -ones(size(mu{2}))};
% 
% fun4Diff_beta = @(mu,beta) [mu{1}.^4, mu{1}.^3, mu{1}.^2, mu{1}];

%% Set the function of parameter constraints

fun = fun1;
funDiff_mu = fun1Diff_mu;
funDiff_beta = fun1Diff_beta;

start = start1;

%% Set the derivatives of the function of parameter constraints (optional)

options.funDiff_mu   = funDiff_mu;
options.funDiff_beta = funDiff_beta;

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
title('Area Calibration Data: Observed vs. fitted values')