function [ParameterMCM,X,A,mu] = funMonteCarloMethod(data,SqrtmU)
% funMonteCarloMethod - Auxiliary function for generating one realization
%  of parameter estimated nbased on using the Monte Carlo Method. 
%
%  For given mean vector 'data' and square root matrix 'SqrtmU' of the
%  specified uncertainty matrix U, the function funMonteCarloMethod generates
%  new data vector mu, which is further used to construct the desing matrix
%  A and the vector of  observations X. Based on that, the model parameters
%  of interest, in particular (MR, M1, M2, Vr, V1, V2), are  estimated as
%  'ParameterMCM' by the Least quares method (LSE). For more details see
%  the Mass Calibration EXAMPLE by Zuda (2023).
%
%  Notation:
%
%  mu(1)  = X1 = X1R_1  % comparison measurement: M1 - MR in evironment 1
%  mu(2)  = X2 = X2R_1  % comparison measurement: M2 - MR in evironment 1
%  mu(3)  = X3 = X21_1  % comparison measurement: M2 - M1 in evironment 1
%  mu(4)  = X4 = X1R_2  % comparison measurement: M1 - MR in evironment 2
%  mu(5)  = X5 = X2R_2  % comparison measurement: M2 - MR in evironment 2
%  mu(6)  = X6 = X21_2  % comparison measurement: M2 - M1 in evironment 2
%  mu(7)  = X7 = X1R_3  % comparison measurement: M1 - MR in evironment 3
%  mu(8)  = X8 = X2R_3  % comparison measurement: M2 - MR in evironment 3
%  mu(9)  = X9 = X21_3  % comparison measurement: M2 - M1 in evironment 3
%  mu(10) = dMR         % difference from nominal of the reference weight
%  mu(11) = VR          % volume of the reference weight
%  mu(12) = G           % gravitational ratio
%  mu(13) = h1 = hR     % center of gravity of the reference weight
%  mu(14) = h2 = h1     % center of gravity of the measured weight 1
%  mu(15) = h3 = h2     % center of gravity of the measured weight 2
%  mu(16) = K           % adjustment constant
%  mu(17) = α           % coefficient of thermal expansion
%  mu(18) = dT1 = dT_1  % temperature differnce T_1 - 20 in evironment 1
%  mu(19) = dT2 = dT_2  % temperature differnce T_2 - 20 in evironment 2
%  mu(20) = dT3 = dT_3  % temperature differnce T_3 - 20 in evironment 3
%  mu(21) = ρ1  = ρ_1   % density of air in evironment 1
%  mu(22) = ρ2  = ρ_2   % density of air in evironment 2
%  mu(23) = ρ3  = ρ_3   % density of air in evironment 3
%
%  Variable related to the vector parameter beta:
%
%  beta(1) = dM1        % difference from nominal of the measured weight 1
%  beta(2) = V1         % volume of the measured weight 1
%  beta(3) = dM2        % difference from nominal of the measured weight 2
%  beta(4) = V2         % volume of the measured weight 2
%
% EXAMPLE
%  data  = [-9.7550e-08 -9.0350e-08  6.3200e-09 -4.7849e-07 9.4810e-08 ...
%            3.8381e-07 -3.2604e-07 -9.3040e-08  2.3331e-07 1+6.4100e-07 ...
%            1.2544e-04  3.0000e-07  3.6500e-02  3.6500e-02 3.0000e-02 ...
%            9.9986e-01  4.8000e-05  8.0000e-01  4.0000e-01 9.0000e-01 ...
%            1.1552e+00  5.8770e-01  8.1380e-01]';
%  ux    = [ 5.0000e-10  3.0000e-10  2.0000e-10  3.0000e-10 3.0000e-10 ...  
%            3.0000e-10  2.0000e-10  3.0000e-10  2.0000e-10 4.0000e-08 ...
%            1.0000e-09  3.0000e-08  5.0000e-04  5.0000e-04 5.0000e-04 ... 
%            5.0000e-06  2.0000e-06  1.0000e-01  1.0000e-01 1.0000e-01 ...
%            5.0000e-04  5.0000e-04  5.0000e-04]';
%  U      = diag(ux.^2);
%  SqrtmU = sqrtm(U)'; 
%  rng(0)
%  N = 1000000;
%  mu = zeros(23,N);
%  ParameterMCM = zeros(6,N);
%  for i = 1:N
%     ParameterMCM(:,i) = funMonteCarloMethod(data,SqrtmU);
%  end
%  BestEstimate = mean(ParameterMCM,2);
%  U_BestEstimate = cov(ParameterMCM');
%  u_BestEstimate = sqrt(diag(U_BestEstimate));
%  alpha = 0.05;
%  coverageFactor = norminv(1-alpha/2);
%  TABLE_MCM = table;
%  TABLE_MCM.Properties.Description = 'Estimated Mass and Volume of the Weights by MCM';
%  TABLE_MCM.ESTIMATE = BestEstimate;
%  TABLE_MCM.STD      = u_BestEstimate;
%  TABLE_MCM.FACTOR   = coverageFactor*ones(size(BestEstimate));
%  TABLE_MCM.LOWER    = BestEstimates - coverageFactor*u_BestEstimate;
%  TABLE_MCM.UPPER    = BestEstimates + coverageFactor*u_BestEstimate;
%  TABLE_MCM.Properties.RowNames = {'M_R (51699)' 'M_1 (15N)' 'M_2 (15963)' 'V_R (51699)' 'V_1 (15N)' 'V_2 (15963)'};
%  TABLE_MCM.Properties.VariableNames = {'ESTIMATE [kg | m^3]'  'STD [kg | m^3]'  'FACTOR'  'LOWER BOUND [kg | m^3]'  'UPPER BOUND [kg | m^3]' };
%  disp(TABLE_MCM)

% Viktor Witkovsky (witkovsky@savba.sk)
% Ver.: 17-Feb-2024 17:02:15

%% ALGORITHM

% Generate new data by the Monte Carlo Method (MCM)
% (23-dimensional vector mu specified my the best estimates - data, 
% and their covariance matrix U

mu = data(:) + SqrtmU * randn(23,1);

% New observation vector (see equation (8))
X = zeros(11,1);
X(1:9) = mu(1:9) * mu(16);
X(10)  = mu(10);
X(11)  = mu(11);

% Realization of random matrix A  (created from equation (8))
A = [ ...
    -(1-mu(12)*mu(13))  (1-mu(12)*mu(14))                 0   (1+mu(17)*mu(18))*mu(21)  -(1+mu(17)*mu(18))*mu(21)                         0; ...
    -(1-mu(12)*mu(13))                  0 (1-mu(12)*mu(15))   (1+mu(17)*mu(18))*mu(21)                          0  -(1+mu(17)*mu(18))*mu(21); ...
                     0 -(1-mu(12)*mu(14)) (1-mu(12)*mu(15))                          0   (1+mu(17)*mu(18))*mu(21)  -(1+mu(17)*mu(18))*mu(21); ...
    -(1-mu(12)*mu(13))  (1-mu(12)*mu(14))                 0   (1+mu(17)*mu(19))*mu(22)  -(1+mu(17)*mu(19))*mu(22)                         0; ...
    -(1-mu(12)*mu(13))                  0 (1-mu(12)*mu(15))   (1+mu(17)*mu(19))*mu(22)                          0  -(1+mu(17)*mu(19))*mu(22); ...
                     0 -(1-mu(12)*mu(14)) (1-mu(12)*mu(15))                          0   (1+mu(17)*mu(19))*mu(22)  -(1+mu(17)*mu(19))*mu(22); ...                    
    -(1-mu(12)*mu(13))  (1-mu(12)*mu(14))                 0   (1+mu(17)*mu(20))*mu(23)  -(1+mu(17)*mu(20))*mu(23)                         0; ...
    -(1-mu(12)*mu(13))                  0 (1-mu(12)*mu(15))   (1+mu(17)*mu(20))*mu(23)                          0  -(1+mu(17)*mu(20))*mu(23); ...
                     0 -(1-mu(12)*mu(14)) (1-mu(12)*mu(15))                          0   (1+mu(17)*mu(20))*mu(23)  -(1+mu(17)*mu(20))*mu(23); ...
                     1                  0                 0                          0                          0                         0; ...
                     0                  0                 0                          1                          0                         0];

ParameterMCM = A \ X;

end