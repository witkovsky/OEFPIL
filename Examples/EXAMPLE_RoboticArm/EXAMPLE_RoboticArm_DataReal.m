%% EXAMPLE_RobotArmSetting_DataReal
%
% This example is inspired by the paper [1], submitted for publication in
% Measurement Science Review by Krepelka, Schovanek, Tucek, Hrabovsky, and
% Jane in 2023 (accepted for publication in 2024).
%
% Here we use the real 3D data presented in [1] (see Table 3).

% Viktor Witkovsky (witkovsky@savba.sk)
% Ver.: 07-Jan-2024 12:29:56

clear 
close all

%% Load Ellipse data

load DataReal.mat

%% Ellipse Data - Table 2 from the paper Krepelka etal MSR 2024

% Measured data 
XYZdata = [ ...
   3514.52                   -579.57                    687.99
   3483.87                   -626.31                    684.07
   3444.28                   -667.52                    678.19
   3397.32                   -701.62                    670.46
   3340.54                   -731.09                    660.67
   3289.08                    -751.1                    651.66
   3275.88                    -763.8                    628.36
    3301.2                   -772.43                    578.81
   3326.56                   -777.42                    527.37
   3350.02                   -779.08                    478.43
   3360.75                   -777.46                    467.53
   3414.21                   -751.57                    467.42
    3462.5                   -716.63                    467.42
   3505.29                      -672                     467.4
   3536.87                   -625.08                    467.33
   3562.01                   -574.14                    467.44
   3583.49                   -516.08                    467.57
   3583.05                   -505.49                    472.19
   3575.69                   -519.99                    522.13
   3561.52                   -536.88                    581.77
   3544.12                   -550.97                    633.32
   3520.35                   -564.97                    685.18];

% Target position
Target = [ ...
   3514.45                   -579.54                    688.22
   3483.84                    -626.3                     684.2
   3444.28                   -667.51                    678.22
    3397.3                   -701.61                    670.52
   3340.52                   -731.08                    660.75
   3289.08                    -751.1                    651.64
   3276.04                   -763.79                    628.44
   3301.33                   -772.42                    578.88
   3326.72                   -777.41                    527.45
   3350.23                   -779.07                    478.53
   3360.75                   -777.47                    467.75
   3414.22                   -751.58                    467.71
   3462.51                   -716.64                     467.7
   3505.31                   -672.01                     467.7
   3536.89                   -625.09                     467.7
   3562.02                   -574.15                     467.7
    3583.5                   -516.08                     467.7
   3583.05                   -505.45                     472.2
   3575.69                   -519.87                    522.16
   3561.53                   -536.74                    581.81
   3544.12                   -550.99                    633.32
   3520.35                   -565.06                    685.16];

[m,n] = size(XYZdata);

%% Fit the data by OEFPIL

% Set the uncertainty matrix.
% Since we have no specific knowledge about the uncertainty structure, the
% uncertainty matrix of measurements U is assumed to be very simple as U =
% sigma2 * I, and the scalar variance sigma2 is subsequently estimated from
% residuals. 
% If U is set as an empty matrix, the algorithm uses the identity matrix.
% Additionally, if options.isEstimatedVariance is set to true, the
% uncertainty matrix is assumed to be U = sigma2 * I, and the scalar
% variance sigma2 is subsequently estimated from residuals.
U = []; 

% Set the required vector of constraints 
% using pre-defined function funHELLA
fun = @(mu,beta) funHELLA(mu,beta,Target);

mu0 = [];                % starting value of the parameters in mu, if empty  
                         % the algorithm uses observed data as mu0
beta0  = [0 0 0 0 0 0]'; % starting value of the parameters in beta

options.q = m*n;          % set the number of restrictions q 

options.tol = 1e-12;     % set the parameter tol for higher precision

options.isEstimatedVariance = 1; % as we have no specific knowledge about 
                                 % the uncertainty matrix of measurement U
                                 % we assume U = sugma2* I, and the scalar 
                                 % variance sigma2 is subsequently
                                 % estimated from residuals.    
options.criterion = 'function';  % set the criterion for optimization 
                                 % (minimize the function value) 

% Fit the data by OEFPIL algorithm                                 
result = OEFPIL(XYZdata,U,fun,mu0,beta0,options);

%% Estimated parameters
beta = result.beta;   % Estimated rotation angles and translation vector
mu = result.mu;       % Estimated mean values of the observed data


% Estimated translation vector [mm]
translationEstimated = beta(4:6);
disp(table(translationEstimated))

% Estimated angles [RAD]
anglesEstimatedRad = beta(1:3);
disp(table(anglesEstimatedRad))

% Estimated angles [DEG]
anglesEstimatedDeg = rad2deg(anglesEstimatedRad);
disp(table(anglesEstimatedDeg))

residuals = result.residuals;
TargetResiduals = reshape(fun(mu,beta),m,n);
EstimatedTarget = Target + TargetResiduals;

%% Compare the transformed data with the target values
% 
% In this context, we assume that the observed data were generated
% through a rotation and translation process from the original target
% values. The mean values of the measurements, mu = [mu_x mu_y mu_z], are
% given by the equation:
%
%   mu = Target * R_xyz' + ones(m, 1) * translation'
%
% where 'translation' is a 3-dimensional column vector, R_xyz = R_z*R_y*R_x 
% is the rotation matrix used, with the first rotation (alpha) about the
% x-axis, the second rotation (beta) about the y-axis, and the third
% rotation (gamma) about the z-axis.
%
% It's important to note that the rotation of a 3D column vector, xyzRot,
% is given by Rxyz * xyz.
%
% However, in our case, we deal with transformed equations because the data
% and mu are represented as (m x 3) matrices, with the 3D coordinates of
% the measured points provided in rows. Consequently, the reverse
% transformation to obtain the Target values from the means of observed
% data involves the following steps:
%
%   (mu - ones(m, 1) * translation') = Target * R_xyz'
%   (mu - ones(m, 1) * translation') = Target * (R_z * R_y * R_x)'
%   (mu - ones(m, 1) * translation') = Target * R_x' * R_y' * R_z'
%
% To identify the Target, we multiply both sides of the equation by the
% inverse rotations. Notably, R_x' * R_x = eye(3), R_y' * R_y = eye(3), and
% R_z' * R_z = eye(3). Thus, we obtain:
%
%   (mu - ones(m, 1) * translation') * R_z * R_y * R_x 
%                                    = Target * R_x' * R_y' * R_z'
%   (mu - ones(m, 1) * translation') * R_z * R_y * R_x = Target
% or
%   Target = (mu - ones(m, 1) * translation') * Rz_ * R_y * R_x
%   Target = (mu - ones(m, 1) * translation') * R_xyz

beta = result.beta;   % Estimated rotation angles and translation vector
mu = result.mu;       % Estimated mean values of the observed data
m = result.m;

a = beta(1);
b = beta(2);
g = beta(3);

ca = cos(a);
cb = cos(b);
cg = cos(g);

sa = sin(a);
sb = sin(b);
sg = sin(g);

R_x = [1  0   0; ...
      0  ca sa; ...
      0 -sa ca];

R_y = [cb  0  -sb; ...
      0   1    0; ...
      sb  0   cb];

R_z = [cg  sg  0; ...
      -sg cg  0; ...
      0    0  1];

R_xyz = R_z*R_y*R_x; % rotation matrix as defined in funHella
translation = beta(4:6);

RotationMatrix = R_xyz;
Translation = ones(m,1)*translation';

%%  Compare XYZmeanTransformed coordinate points with the target points

Means_XYZ  = mu;
TransformedMeans = (Means_XYZ - Translation) * RotationMatrix;

Tab1 = table(Target);
Tab2 = table(TransformedMeans);

disp(Tab1) 
disp(Tab2)

TransformedMeans_Residuals = TransformedMeans - Target;

Tab3 = table(TransformedMeans_Residuals);
disp(Tab3)

%%  Compare TransformedData coordinate points with the target points

Data_XYZ  = XYZdata;
TransformedData = (Data_XYZ - Translation) * RotationMatrix;

Tab1 = table(Target);
Tab4 = table(TransformedData);
disp(Tab1) 
disp(Tab4)


TransformedData_Residuals = TransformedData - Target;
Tab5 = table(TransformedData_Residuals);
disp(Tab5)

%%  Compare TransformedData coordinate points evaluated from OEFPIL fit 
%   with the points calculate by authors in [1]  

% In the paper [1], the authors denoted their transformed points as
% **corrected** points and in Table 3 are presented the following results:
% Corrected = [ ...
%         3514.49                  -579.507                   688.114
%         3483.853                  -626.256                   684.198
%         3444.275                  -667.478                   678.317
%         3397.326                  -701.592                   670.582
%         3340.557                  -731.079                   660.783
%         3289.105                  -751.105                   651.762
%         3275.915                  -763.815                   628.462
%         3301.253                  -772.452                   578.922
%         3326.63                  -777.451                   527.492
%         3350.106                  -779.119                   478.559
%         3360.839                  -777.499                   467.662
%         3414.293                  -751.596                   467.562
%         3462.574                  -716.643                   467.567
%         3505.352                  -672.003                   467.548
%         3536.921                  -625.075                   467.474
%         3562.048                  -574.128                   467.577
%         3583.513                  -516.063                   467.698
%         3583.069                  -505.472                   472.314
%         3575.697                  -519.959                   522.256
%         3561.512                  -536.836                   581.897
%         3544.1                    -550.915                   633.445
%         3520.317                  -564.907                   685.302];

% For comparison, see the transformed OEFPIL data (rounded to three decimal
% places) given in TransformedData 
% TransformedData = [ ...
%         3514.488                  -579.508                   688.113
%         3483.852                  -626.257                   684.197
%         3444.274                  -667.479                   678.316
%         3397.325                  -701.593                    670.58
%         3340.556                   -731.08                    660.78
%         3289.104                  -751.105                    651.76
%         3275.914                  -763.815                   628.459
%         3301.252                  -772.453                    578.92
%          3326.63                  -777.451                   527.489
%         3350.106                  -779.119                   478.557
%         3360.839                  -777.499                    467.66
%         3414.293                  -751.596                    467.56
%         3462.574                  -716.644                   467.565
%         3505.352                  -672.003                   467.546
%         3536.921                  -625.075                   467.473
%         3562.048                  -574.129                   467.577
%         3583.513                  -516.063                   467.697
%         3583.069                  -505.472                   472.314
%         3575.696                   -519.96                   522.256
%         3561.512                  -536.836                   581.896
%         3544.099                  -550.916                   633.444
%         3520.316                  -564.907                   685.301];

% Here is the difference between the OEFPIL and original solution in [1]:

Differences = TransformedData - Corrected;
Tab6 = table(Differences);
disp(Tab6)
