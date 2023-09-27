%% EXAMPLE: Flow Meter / Fractional Polynomial Calibration

clear
close all

%% Data

x   = [3.1100e+03 2.5000e+03 1.5000e+03 7.4000e+02 4.2200e+02 ...
    1.2000e+02 2.5000e+01 1.0000e+01 7.0000e+00 5.5000e+00 ...
    3.2500e+00]';

y   = [-1.4040e+00 -1.1500e+00 -6.0900e-01 2.6200e-01 9.8000e-01 ...
    1.9290e+00 1.3400e+00 -2.5400e-01 -1.3800e+00 -2.4220e+00 ...
    -5.2990e+00]';

%% Uncertainty matrix

m   = 11;

Ux1 = [9.7000e+00 3.1000e+00 1.1000e+00 2.7000e-01 8.7000e-02 7.1000e-03;
    3.1000e+00 6.3000e+00 1.1000e+00 2.7000e-01 8.7000e-02 7.1000e-03;
    1.1000e+00 1.1000e+00 2.3000e+00 2.7000e-01 8.7000e-02 7.1000e-03;
    2.7000e-01 2.7000e-01 2.7000e-01 5.5000e-01 8.7000e-02 7.1000e-03;
    8.7000e-02 8.7000e-02 8.7000e-02 8.7000e-02 1.8000e-01 7.1000e-03;
    7.1000e-03 7.1000e-03 7.1000e-03 7.1000e-03 7.1000e-03 1.4000e-02];

Ux2 = [2.5000e-03 2.0000e-04 9.6000e-05 5.9000e-05 2.1000e-05;
    2.0000e-04 4.0000e-04 9.6000e-05 5.9000e-05 2.1000e-05;
    9.6000e-05 9.6000e-05 2.0000e-04 5.9000e-05 2.1000e-05;
    5.9000e-05 5.9000e-05 5.9000e-05 1.2000e-04 2.1000e-05;
    2.1000e-05 2.1000e-05 2.1000e-05 2.1000e-05 4.2000e-05];
Ux  = blkdiag(Ux1,Ux2);
Uy  = diag([1.0240e-03 3.7210e-03 2.6010e-03 1.7640e-03 3.7210e-03 ...
    4.0960e-03 3.3640e-03 3.2490e-03 6.0840e-03 6.2410e-03 ...
    1.0609e-02]);
Uxy = zeros(m);

%% Function - nonlinear constraints on the parameter

fun  = @(mu,nu,beta) beta(1)*mu.^(-1) + beta(2)*mu.^(-0.5) + ...
    beta(3) + beta(4)*mu.^(0.5) + beta(5)*mu - nu;

%% Staring values

mu0   = x;
nu0   = y;
beta0 = [0; 0; 0; 0; 1];

%% Optional settings / method and criterion function

options.method = 'oefpil';
options.criterion = 'parameterdifferences';
options.tol = 1e-10;

%% Calculate the result

result = OEFPIL2D(x,y,{Ux,Uy,Uxy},fun,mu0,nu0,beta0,options);

disp(result)