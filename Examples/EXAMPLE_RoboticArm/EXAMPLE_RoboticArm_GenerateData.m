%% EXAMPLE_RobotArmSetting_GenerateData
%
% Generate and fit artificial data inspired by the paper [1], submitted for
% publication in Measurement Science Review by Krepelka, Schovanek, Tucek,
% Hrabovsky, and Jane in 2023 (accepted for publication in 2024).
% 
% For details see also, EXAMPLE_RoboticArm_OEFPIL.m
%
% REFERENCES
% [1]  Krepelka J., Schovanek P., Tucek P., Hrabovsky M., Jane F.
%      Optimization of component assembly in automotive industry.
%      Measurement Science Review 2024.

% Viktor Witkovsky (witkovsky@savba.sk)
% Ver.: 07-Jan-2024 12:29:56

clear 
close all
%% Generate artificial data

% Set the true rotation angles alpha, beta, gamma in radians 
% If the angles are set in degrees use transformation to radians 
% beta1 = deg2rad([-8 -5 -15]);
angle = [-0.1396   -0.0873   -0.2618]';

% Set the true translations in direction x, y, and z in mm
translation = [-1 -15 -20]';

% For simplicity specify the angles a (alpha), b (beta), and g (gamma)
a = angle(1);
b = angle(2);
g = angle(3);

% Pre-calculate the cosine and sine values
ca = cos(a);
cb = cos(b);
cg = cos(g);

sa = sin(a);
sb = sin(b);
sg = sin(g);

% Rotation matrix R_zyx =  R_x*R_y*R_z (first rotation with angle g (gamma)
% about z-axis, second rotation with angle b (beta) about y-axis, third
% rotation with angle a (alpha) about x-axis)  

R_zyx = [           cb*cg,            cb*sg,   -sb; ...
         cg*sa*sb - ca*sg, ca*cg + sa*sb*sg, cb*sa; ...
         sa*sg + ca*cg*sb, ca*sb*sg - cg*sa, ca*cb];

% The reverse rotation matrix R_xyz =  R_z*R_y*R_x (first rotation with
% angle a (alpha) about x-axis, second rotation with angle b (beta) about
% y-axis, third rotation with angle g (gamma) about z-axis) 
% (Note the error in R_xyz(2,3) in the paper [1])
% R_xyz  = [cb*cg ca*sg+sa*sb*cg sa*sg-ca*sb*cg; ...
%          -cb*sg ca*cg-sa*sb*sg sa*sg+ca*sb*sg; ...
%           sb -sa*cb ca*cb];

R_xyz = [ cb*cg, ca*sg + cg*sa*sb, sa*sg - ca*cg*sb; ...
         -cb*sg, ca*cg - sa*sb*sg, cg*sa + ca*sb*sg; ...
             sb,           -cb*sa,            ca*cb];

R_x = [1  0   0; ...
       0  ca sa; ...
       0 -sa ca];

R_y = [cb  0  -sb; ...
       0   1    0; ...
       sb  0   cb];

R_z = [cg  sg  0; ...
       -sg cg  0; ...
       0    0  1];

Target = [ ...
   590    10   175
   560    10   175
   590   -10   175
   560   -10   175
   590    25   140
   560    25   140
   590    25   160
   560    25   160
   600   -15   135
   600    15   135
   600   -15   165
   600    15   165];

[m,n] = size(Target);

% Expectations mu of the measurements derived by rotating and translating
% the 'target' values Target
% mu = Target * (R_z*R_y*R_x)' + ones(m,1)*translation'; 

mu = Target*R_xyz' + ones(m,1)*translation'; 

% Check correctness. Backward transformation - Target values received by
% the transformation from the mean values of the measurements
%Target = (mu - ones(m,1)*translation')* (Rx'*Ry'*Rz')';

Target = (mu - ones(m,1)*translation')* R_xyz;

%f2 = (mu - ones(m,1)*beta2')* Rzyx' - XYZ0;

%% Generate data

rng(0) % for reproducibility

% Set the true values of standard deviations sigma of the block diagonal
% matrices 
sigma = [1 1 1]/100;

X = mu(:,1) + sigma(1)*randn(m,1);
Y = mu(:,2) + sigma(2)*randn(m,1);
Z = mu(:,3) + sigma(3)*randn(m,1);

XYZ = [X Y Z];

% XYZ = [ ...]
%   585.1827  127.1157  102.5871
%   556.3282  119.3728  105.2082
%   590.0464  107.9222   99.8171
%   561.2101  100.1780  102.4225
%   577.3323  145.4243   70.1368
%   548.4486  137.7054   72.7417
%   579.7130  143.1984   89.8728
%   550.8534  135.4635   92.4681
%   596.1736  110.1841   58.7741
%   588.8280  138.9554   62.9352
%   599.7067  106.8228   88.3504
%   592.4131  135.6220   92.5527];


%% Fit the EIV model by OEFPIL

% Set the observed data
data = XYZ;

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

% Alternatively, if we have qualified information about the uncertaity
% matrix U, we can specify it as a cell structure (for details see OEFPIL) 
% sigma = [1 1 1]/100;
% sigma2 = sigma.^2;
% U = { sigma2(1)*eye(m), sigma2(2)*eye(m), sigma2(3)*eye(m) }
   
% Set the required vector of constraints 
% using pre-defined function funHELLA
fun = @(mu,beta) funHELLA(mu,beta,Target);

mu0 = [];                % starting value of the parameters in mu, if empty  
                         % the algorithm uses observed data as mu0
beta0  = [0 0 0 0 0 0]'; % starting value of the parameters in beta

options.q = 36;          % set the number of restrictions q 

options.tol = 1e-12;     % set the parameter tol for higher precision

options.isEstimatedVariance = 1; % as we have no specific knowledge about 
                                 % the uncertainty matrix of measurement U
                                 % we assume U = sugma2* I, and the scalar 
                                 % variance sigma2 is subsequently
                                 % estimated from residuals.    
options.criterion = 'function';  % set the criterion for optimization 
                                 % (minimize the function value) 

% Fit the data by OEFPIL algorithm                                 
result = OEFPIL(data,U,fun,mu0,beta0,options);

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

%% Compare result 

beta = result.beta;   % Estimated rotation angles and translation vector
mu   = result.mu;     % Estimated mean values of the observed data
m    = result.m;

% For simplicity specify the estimated angles a (alpha), b (beta), and g
% (gamma)
a = beta(1);
b = beta(2);
g = beta(3);

% Pre-calculate the cosine and sine values
ca = cos(a);
cb = cos(b);
cg = cos(g);

sa = sin(a);
sb = sin(b);
sg = sin(g);

% Specify the rotation matrices with estimated angles
R_x = [1  0   0; ...
       0  ca sa; ...
       0 -sa ca];

R_y = [cb  0  -sb; ...
       0   1    0; ...
       sb  0   cb];

R_z = [cg  sg  0; ...
       -sg cg  0; ...
       0    0  1];

R_xyz = R_z*R_y*R_x; % rotation matrix with estimated angles 
                     % as defined in funHella

% Specify the estimated translation vector
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

Data_XYZ  = data;
TransformedData = (Data_XYZ - Translation) * RotationMatrix;

Tab4 = table(TransformedData);
disp(Tab1) 
disp(Tab4)


TransformedData_Residuals = TransformedData - Target;
Tab5 = table(TransformedData_Residuals);
disp(Tab5)