%% EXAMPLE_RobotArmSetting_OEFPIL
%
% This example is inspired by the paper [1], submitted for publication in
% Measurement Science Review by Krepelka, Schovanek, Tucek, Hrabovsky, and
% Jane in 2023 (accepted for publication in 2024):
%
% - The paper addresses the positioning of glued parts by robots in the
%   manufacturing process of automotive headlights, requiring high precision
%   and efficiency.
% - The authors propose a mathematical method to optimize the robot arm
%   setting, ensuring the closest contact of the glued parts based on the
%   position of the nominal and control points on the contact surfaces of the
%   two joined parts.
% - The method involves finding the minimum of an objective function that
%   combines the sum of squared deviations of the position of the control
%   points from the nominal points and the sum of squared deviations from the
%   root mean square deviation.
% - The method can be applied to different types of transformations, such as
%   translation followed by rotation or rotation followed by translation, and
%   different orders of rotation axes.
% - The paper presents numerical simulations and a practical example to
%   demonstrate the effectiveness and accuracy of the proposed method.
%
% Motivated by the paper [1], the example RobotArmSetting_OEFPIL provides
% an illustration of estimating parameters of interest and their
% uncertainties for a robot arm setting. The objective is to ensure the
% contact of assembly parts based on the initial control position (measured
% data) and the nominal points.
%
% In this scenario, the data consists of measured 3D coordinates (along
% with a specified uncertainty matrix) of points on the assembly part
% (e.g., automotive headlight part) at its initial (starting) position and
% the given target coordinates.
%
% The assumption is that the robot arm can reach the target by making
% translation followed by rotation (or rotation followed by translation).
% The required translation is specified by a 3-dimensional vector, and the
% required rotations (about axes x, y, and z) are specified by the rotation
% angles alpha, beta, and gamma. Thus, the parameters of interest are the
% three rotation angles and the 3-dimensional translation vector.
%
% For solving this problem, the example suggests using a regression model of
% direct measurements with nonlinear constraints on its parameters. It is an
% Errors-in-Variable (EIV) model, as the required relations are expressed by
% the nonlinear function of the mean values of the measurements and the
% model parameters of interest.
%
% The positions of the initial coordinates of the specified m points on the
% assembly part are modeled as direct measurements (with errors): 
%   X = mu_x + epsilon_x, Y = mu_y + epsilon_y, Z = mu_z + epsilon_z.
% The complete knowledge about uncertainties is expressed by the uncertainty
% (3*m x 3*m)-dimensional matrix U, with the structure specified for each
% block of the matrix.
%
% Important parts of the model are the required (nonlinear) constraints on
% its unknown parameters mu_x, mu_y, mu_z, alpha, beta, gamma, and the
% vector of translation (the required translations in each direction)
% specified by translation_x, translation_y, and translation_z. 
% 
% The following formulation requires that the expected values of the
% measured coordinates are first translated and then rotated about z-axis,
% y-axis, and x-axes (in this specific ordering), such that such
% transformed coordinates are equal to the targeted values.
%
% In fact, the q-dimensional vector of restrictions, f, is specified as:
%   f = vec( ([mu_x, mu_y, mu_z] 
%       - ones(m, 1) * [translation_x, translation_y, translation_z]') * R_zyx  
%       - [target_x, target_y, target_z] )
% where R_zyx = R_x * R_y * R_z is the rotation matrix with the first
% rotation with angle gamma about the axis z, the second rotation with
% angle beta about the axis y, and the third rotation with angle alpha
% about the axis x). The (m x 3)-matrix XYZ0 = [target_x, target_y,
% target_z] is a given matrix of target coordinates.
%
% Here, the OEFPIL algorithm [2] is presented as a suitable tool for fitting
% measured data and estimating the required parameters of interest together
% with their uncertainties.
%
% REFERENCES
% [1]  Krepelka J., Schovanek P., Tucek P., Hrabovsky M., Jane F.
%      Optimization of component assembly in the automotive industry.
%      Measurement Science Review 2024.
% [2]  Charvatova Campbell, A., Slesinger, R., Klapetek, P., Chvostekova,
%      M., Hajzokova, L., Witkovsky, V. and Wimmer, G. Locally best linear
%      unbiased estimation of regression curves specified by nonlinear
%      constraints on the model parameters. AMCTMT 2023 - Advanced
%      Mathematical and Computational Tools in Metrology and Testing 2023
%      Sarajevo, Bosnia and Herzegovina, 26-28 September 2023.
%

% Viktor Witkovsky (witkovsky@savba.sk)
% Ver.: 07-Jan-2024 12:29:56

clear 
close all

%% Data from the paper Krepelka etal MSR 2024

% First measured data /Intermediate position of the robotic arm P1
data_Position1 = [ ...
  591.5250    6.8350  172.5430
  561.6550    9.4660  173.4510
  589.8180  -13.0400  173.9770
  559.9480  -10.4100  174.8850
  591.5330   19.3330  136.5730
  561.6620   21.9640  137.4810
  592.2600   20.7090  156.5130
  562.3890   23.3400  157.4210
  597.8940  -21.6390  134.1540
  600.4540    8.1740  132.0030
  598.9850  -19.5740  164.0630
  601.5450   10.2390  161.9120];

% Second measured data / Intermediate position of the robotic arm P2
data_Position2 = [ ...
  590.5300    8.7260  174.2960
  560.5490    9.7730  174.5310
  589.8360  -11.2550  174.8260
  559.8550  -10.2080  175.0600
  590.7440   22.7960  138.9130
  560.7640   23.8420  139.1470
  590.9190   23.3190  158.9050
  560.9380   24.3660  159.1390
  599.3080  -17.6460  134.8950
  600.3480   12.3250  134.1010
  599.5690  -16.8610  164.8840
  600.6090   13.1110  164.0900];

% Final position / Target position of the robotic arm 
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

%% Fit the data by OEFPIL

data = data_Position1;

% Alternatively choose the second measured data 
% Intermediate position of the robotic arm P2
% data = data_Position2;

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

Data_XYZ  = data;
TransformedData = (Data_XYZ - Translation) * RotationMatrix;

Tab4 = table(TransformedData);
disp(Tab1) 
disp(Tab4)


TransformedData_Residuals = TransformedData - Target;
Tab5 = table(TransformedData_Residuals);
disp(Tab5)