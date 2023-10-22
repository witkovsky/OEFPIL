%% EXAMPLE_EllipseFit_Data
%  Generate artificial data

clear
close all
 
%% Data

alpha0 = 0; beta0 = 0;        % true ellipse center [0,0]
alpha1 = 1; beta1 = 0.75;     % true amplitudes
phi0 = pi/3;                  % phase shift
X = @(t) alpha0 + alpha1 * cos(t);
Y = @(t) beta0 + beta1 * sin(t + phi0);

m       = 15;                 % No. of observations (x,y)
Ncycles = 0.8;                % No. of whole ellipse cycles
sigma   = 0.05;               % true error STD
phi     = Ncycles*(2*pi)*sort(rand(m,1)); % true phases

mu = X(phi);
nu = Y(phi);

% munu = [ ...
%     0.9970    0.6767
%     0.9651    0.7250
%     0.8608    0.7500
%     0.6673    0.7127
%     0.5438    0.6679
%     0.2163    0.5066
%    -0.1353    0.2837
%    -0.2966    0.1655
%    -0.5266   -0.0233
%    -0.7775   -0.2691
%    -0.9258   -0.4595
%    -0.8434   -0.7493
%    -0.0862   -0.4296
%     0.0710   -0.3280
%     0.2319   -0.2141];

x = mu + sigma*randn(m,1);
y = nu + sigma*randn(m,1);

% xy = [ ...
%     0.9923    0.6329
%     0.9820    0.7410
%     0.8156    0.7220
%     0.6528    0.6971
%     0.5613    0.6394
%     0.1245    0.4553
%    -0.0835    0.2383
%    -0.1754    0.1550
%    -0.4786   -0.1082
%    -0.7933   -0.2388
%    -0.9044   -0.4654
%    -0.8952   -0.7143
%     0.0077   -0.4161
%     0.1180   -0.3032
%     0.2713   -0.2883];
