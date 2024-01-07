function [f, fmat, beta, mu, XYZ0, Rxyz, Rzyx, Rx, Ry, Rz] = funHELLA(mu, beta, XYZ0)
% FUNHELLA Computes the q-dimensional vector function for the
% rotation-translation function transformation. This function assumes that
% mu = {mu1, mu2, mu3} is a cell with m-dimensional vectors of mean vector
% parameters in n = 3 dimensions, and beta is a p-dimensional vector of
% parameters of interest. Here, p = 6, and beta = [beta1; beta2], where
% beta1 = [alpha, beta, gamma]' are the rotation angles about the axes x,
% y, and z, and beta2 is a 3-dimensional vector of translation, beta2 =
% [x0, y0, z0]'.
% 
% For more details, refer to the paper "Optimization of component assembly
% in the automotive industry" by Krepelka et al., MSR (2024).
% 
% The q-dimensional vector of restrictions, f, is specified as:
%   f = vec( ([mu1, mu2, mu3] - ones(m, 1) * beta2') * Rzyx  - XYZ0)
% where Rzyx = Rx * Ry * Rz is the rotation matrix with the first rotation
% with angle gamma about the axis z, the second rotation with angle beta
% about the axis y, and the third rotation with angle alpha about the axis
% x).
%
% The (m x 3)-matrix XYZ0 is a given matrix of target coordinates of the 3D
% points after the rotation-translation transformation and is specified as
% XYZ0 = [X0, Y0, Z0].
%
% SYNTAX:
%    [f, fmat, beta, mu, XYZ0, Rxyz, Rzyx, Rx, Ry, Rz] = funHELLA(mu, beta, XYZ0)
%
% INPUTS
%  mu    - is an (m x 3) matrix mu = [mu1, mu2, mu3], or a cell with three
%          m-dimensional vectors, the mean vector parameters, mu =
%          {mu1, mu2, mu3}.
%  beta  - is a p-dimensional vector of parameters of interest, where p = 6,
%          and beta = [beta1; beta2], where beta1 = [alpha, beta, gamma]'
%          are the rotation angles about the axes x, y, and y,
%          respectively, and beta2 is a 3-dimensional vector of
%          translation, beta2 = [x0, y0, z0]'.
%  XYZ0  - is an (m x 3) matrix XYZ0 = [X0, Y0, Z0], or a cell with three
%          m-dimensional vectors, the specified target coordinates X0, Y0,
%          and Z0, namely, XYZ0 = {X0, Y0, Z0}.
%
% Notes: In this context, we assume that the observed data were generated
% through a rotation and translation process from the original target
% values. The mean values of the measurements, mu = [mu1 mu2 mu3], are
% given by the equation:
%
%   mu = Target * Rxyz' + ones(m, 1) * translation'
%
% where 'translation' is a 3-dimensional column vector, Rxyz = Rx * Ry * Rz
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
%   (mu - ones(m, 1) * translation') = Target * Rxyz'
%   (mu - ones(m, 1) * translation') = Target * (Rz * Ry * Rx)'
%   (mu - ones(m, 1) * translation') = Target * Rx' * Ry' * Rz'
%
% To identify the Target, we multiply both sides of the equation by the
% inverse rotations. Notably, Rx' * Rx = eye(3), Ry' * Ry = eye(3), and Rz'
% * Rz = eye(3). Thus, we obtain:
%
%   (mu - ones(m, 1) * translation') * Rz * Ry * Rx = Target * Rx' * Ry' * Rz'
%   (mu - ones(m, 1) * translation') * Rz * Ry * Rx = Target
% or
%   Target = (mu - ones(m, 1) * translation') * Rz * Ry * Rx
%   Target = (mu - ones(m, 1) * translation') * Rxyz
%
% EXAMPLE (Component assembly in automotive industry / Fit by OEFPIL)
% % Measured 3D coordinates of 12 points of the component moved by the
% % robotic arm 
%   XYZ = [ ...
%    585.1773  127.1079  102.5819
%    556.3108  119.3742  105.1960
%    590.0680  107.9160   99.8094
%    561.2000  100.1797  102.4237
%    577.3278  145.4257   70.1346
%    548.4621  137.6893   72.7488
%    579.7173  143.1833   89.8640
%    550.8497  135.4494   92.4790
%    596.1393  110.1782   58.7851
%    588.8004  138.9699   62.9427
%    599.7204  106.8150   88.3801
%    592.3842  135.6058   92.5391];
% % Specified uncertainty matrix U. 
% % If empty, the algorithm OEFFPIL will change and set U = I 
%   U = [];
% % Constraints on the model parameters specific for the required
% % translation-rotation transformation
%   fun = @(mu, beta) funHELLA(mu, beta, XYZ0);
% % Starting values
%   mu0 = XYZ;
%   beta0 = [0 0 0 0 0 0]';
% % Specified controls/options for the algorithm OEFPIL
%   options.tol = 1e-13;
%   options.q = 36; % Number of restrictions on parameters specified by fun
% % as we have no specific knowledge about  the covariance matrix of
% % measurements we assume Sigma = sigma2* U, and the scalar variance
% % sigma2 is subsequently estimated from residuals.   
%   options.isEstimatedVariance = true;
% % Run the algorithm OEFPIL
%   result = OEFPIL(XYZ, U, fun, mu0, beta0, options);
%
% REFERENCES
% [1]  Krepelka J., Schovanek P., Tucek P., Hrabovsky M., Jane F.
%      Optimization of component assembly in the automotive industry.
%      Measurement Science Review 2024.
% [2]  Charvatova Campbell, A., Slesinger, R., Klapetek, P., Chvostekova,
%      M., Hajzokova, L., Witkovsky, V. Locally best linear unbiased
%      estimation of regression curves specified by nonlinear constraints on
%      the model parameters. AMCTMT 2023 - Advanced Mathematical and
%      Computational Tools in Metrology and Testing 2023 Sarajevo, Bosnia and
%      Herzegovina, 26-28 September 2023.

% Viktor Witkovsky (witkovsky@savba.sk)
% Ver.: 07-Jan-2024 12:26:21

%% ALGORITHM
if iscell(mu)
    mu = [mu{1} mu{2} mu{3}];
end

if iscell(XYZ0)
    XYZ0 = [XYZ0{1} XYZ0{2} XYZ0{3}];
end

m = size(mu,1);

beta1 = beta(1:3);
beta2 = beta(4:6);

a = beta1(1);
b = beta1(2);
g = beta1(3);

ca = cos(a);
cb = cos(b);
cg = cos(g);

sa = sin(a);
sb = sin(b);
sg = sin(g);

% Rotation matrix Rx - rotation with angle a (alpha) aroud axis x
Rx = [1  0   0; ...
      0  ca sa; ...
      0 -sa ca];

% Rotation matrix Ry - rotation with angle b (beta) aroud axis y
Ry = [cb  0  -sb; ...
      0   1    0; ...
      sb  0   cb];

% Rotation matrix Rz - rotation with angle g (gamma) aroud axis z
Rz = [cg  sg  0; ...
      -sg cg  0; ...
      0    0  1];

% Rotation matrix Rzyx =  Rx*Ry*Rz (first rotation with angle g (gamma)
% about axis z, second rotation with angle b (beta) about axis y, third
% rotation with angle a (alpha) about axis x)  
% Rzyx  = [cb*cg cb*sg -sb; ...
%          sa*sb*cg-ca*sg ca*cg+sa*sb*sg sa*cb; ...
%          sa*sg+ca*sb*cg ca*sb*sg-sa*cg ca*cb];
Rzyx  = [           cb*cg,            cb*sg,   -sb; ...
         cg*sa*sb - ca*sg, ca*cg + sa*sb*sg, cb*sa; ...
         sa*sg + ca*cg*sb, ca*sb*sg - cg*sa, ca*cb];

% The reverse rotation matrix Rxyz =  Rz*Ry*Rx (first rotation with angle a (alpha)
% about axis x, second rotation with angle b (beta) about axis y, third
% rotation with angle g (gamma) about axis z)  is 
% (Note the error in Rxyz(2,3) in the MSR paper)
% Rxyz  = [cb*cg ca*sg+sa*sb*cg sa*sg-ca*sb*cg; ...
%         -cb*sg ca*cg-sa*sb*sg sa*sg+ca*sb*sg; ...
%          sb -sa*cb ca*cb];
Rxyz  = [ cb*cg, ca*sg + cg*sa*sb, sa*sg - ca*cg*sb; ...
         -cb*sg, ca*cg - sa*sb*sg, cg*sa + ca*sb*sg; ...
             sb,           -cb*sa,            ca*cb];

fmat = (mu - ones(m,1)*beta2')*Rxyz - XYZ0; % Rotation matrix Rxyz
%fmat = (mu - ones(m,1)*beta2')*Rzyx - XYZ0; % Rotation matrix Rzyx

f = fmat(:);

end