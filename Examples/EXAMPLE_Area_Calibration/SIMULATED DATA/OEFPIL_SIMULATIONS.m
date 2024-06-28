%% SIMULATIONS

% Viktor Witkovsky (witkovsky@savba.sk)
% Ver.: '18-Sep-2023 14:57:24'

clear
close all

%% Load data / Good data

%load Data14_BV52_1622021alahysitron.mat
%load Data16_BV52_2222021alahysitron.mat
%load Data17_BV52_932021progressive
%load Data18_BV52_19112020.mat
%load Data19_BV52_322021.mat
%load DataBerk52_25032021.mat

%% Bad data

load DataBerk52_25032021.mat; Uy = Uy20;
%load DataBerk227_2452.mat; Uy = Uy20;
%load DataBerk618_2met.mat

%%

q = length(x);
if isempty(Uxy)
    Uxy = zeros(q);
end
U   = [Ux Uxy; Uxy' Uy];

%% Model / %% Generate data
% EXAMPLE_OEFPIL DataBerk227_2452
%     "params": ["a2", "a1", "a12", "a14"],
%     "f": "a2*x^2 + a1*x + a12*x^(1/2) + a14*x^(1/4)",
%     "dfdx": ["2*a2*x + a1 + 1/2*a12*x^(-1/2) + 1/4*a14*x^(-3/4)"],
%     "dfdp": ["x^2", "x", "x^(1/2)", "x^(1/4)"]

gfun = @(mu,beta)beta(1)*mu.^2+beta(2)*mu+beta(3)*mu.^(1/2)+beta(4)*mu.^(1/4);

beta = [ 5   10   500  -1000];

mu = round(x*10)/10;
nu = gfun(mu,beta);

munu = [mu nu];

xy = munu(:) + sqrtm(full(U)) * randn(2*q,1);

x = xy(1:q);
y = xy(q+(1:q));

%% Observed vs. True Values

figure
plot(x,y,'*',mu,nu,'o')
xlabel('x')
ylabel('y')
title('Observed vs. True Values')
legend('Observed values','True values','Location','nw')

%% Fit the model by OEFPIL

mu0 = {x, y};

beta0 = [1;0;0;0];

fun = @(mu,beta)beta(1)*mu{1}.^2+beta(2)*mu{1}+beta(3)*mu{1}.^(1/2)+beta(4)*mu{1}.^(1/4)-mu{2};

options.funDiff_mu = @(mu,beta) {2*beta(1)*mu{1} + beta(2) + mu{1}.^(-1/2)*beta(3)/2 + mu{1}.^(-3/2)*beta(4)/4,  - ones(size(mu{2}))};
options.funDiff_beta = @(mu,beta) [mu{1}.^2, mu{1}, mu{1}.^(1/2), mu{1}.^(1/4)];

options.method = 'oefpil';
options.tol = 1e-10;
options.criterion = 'parameterdifferences';
%options.criterion = 'function';
options.isPlot = false;

result = OEFPIL({x,y},U,fun,mu0,beta0,options);

%% Observed vs. True vs. Fitted Values

figure
plot(x,y,'*',mu,nu,'o',result.muCell{1},result.muCell{2},'+')
xlabel('x')
ylabel('y')
title('Observed vs. True vs. Fitted Values')
legend('Observed values','True values','Fitted Values','Location','nw')
