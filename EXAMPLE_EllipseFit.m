%% EXAMPLE: Ellipse fit the ellipse for generated measurements x and y)

clear
close all

%% Data

x     = [ 1.0310    1.0150    0.9677    0.8924    0.8025 ...
    -0.3599   -0.7536   -0.7702   -1.0107   -0.9757 ...
    -0.8163   -0.6509   -0.3045   -0.0935    0.3155]';

y     = [ 0.7332    0.6956    0.7028    0.7678    0.7955 ...
    0.1984   -0.3282   -0.3222   -0.7162   -0.7561 ...
    -0.6521   -0.7264   -0.6397   -0.4567   -0.2040]';

%% Uncertainty matrix

m     = 15;
Ux    = 0.05^2*eye(m);
Uy    = 0.05^2*eye(m);
Uxy   = zeros(m);
U     = [Ux Uxy; Uxy' Uy];

%% Function - nonlinear constraints on the parameter

fun   = @(mu,nu,beta) mu.^2 + beta(1)*nu.^2 + beta(2)*mu.*nu + ...
    beta(3)*mu + beta(4)*nu + beta(5);

%% Staring values

mu0   = x;
nu0   = y;
beta0 = [1; 0; 0; 0; -0.1];

%% Optional settings /  method

options.method = 'oefpil';

%% Calculate the result

result = OEFPIL2D(x,y,U,fun,mu0,nu0,beta0,options);

disp(result)