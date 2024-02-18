function fun = funMassCalibration(mu,beta)
% funMassCalibration - defines the necessary restrictions on the EIV model
%  parameters for the calibration of weight masses. 
% 
%  funMassCalibration is an anonymous function that takes arguments mu and
%  beta and returns fun, a q-dimensional column vector of constraints.
%  Ideally, these constraints should equal zero. In this context, mu can be
%  specified as either an m-dimensional column vector or as a cell array
%  with one column vector {mu}, while beta is a p-dimensional vector. The
%  output fun from funMassCalibration(mu, beta) represents the q
%  constraints on the model parameters mu and beta. If the model
%  restrictions are satisfied, fun = funMassCalibration(mu, beta) = 0.
%   
%  Variable related to the vector parameter mu:
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
%  beta(2) = dM2        % difference from nominal of the measured weight 2
%  beta(3) = V1         % volume of the measured weight 1
%  beta(4) = V2         % volume of the measured weight 2
%
% EXAMPLE
%  data  = [-9.7550e-08 -9.0350e-08  6.3200e-09 -4.7849e-07 9.4810e-08 ...
%            3.8381e-07 -3.2604e-07 -9.3040e-08  2.3331e-07 6.4100e-07 ...
%            1.2544e-04  3.0000e-07  3.6500e-02  3.6500e-02 3.0000e-02 ...
%            9.9986e-01  4.8000e-05  8.0000e-01  4.0000e-01 9.0000e-01 ...
%            1.1552e+00  5.8770e-01  8.1380e-01]';
%  ux    = [ 5.0000e-10  3.0000e-10  2.0000e-10  3.0000e-10 3.0000e-10 ...  
%            3.0000e-10  2.0000e-10  3.0000e-10  2.0000e-10 4.0000e-08 ...
%            1.0000e-09  3.0000e-08  5.0000e-04  5.0000e-04 5.0000e-04 ... 
%            5.0000e-06  2.0000e-06  1.0000e-01  1.0000e-01 1.0000e-01 ...
%            5.0000e-04  5.0000e-04  5.0000e-04]';
%  U     = diag(ux.^2);
%  fun   = @(mu,beta) funMassCalibration(mu,beta)
%  mu0   = data;
%  beta0 = [0 0 0 0]';
%  clear options
%  options.method        = 'oefpilrs2';
%  options.criterion     = 'function';
%  options.numDiffMethod = 'LynnesMoller';
%  options.q             = 9;
%  options.tol           = 1e-15;
%  result = OEFPIL(data,U,fun,mu0,beta0,options);

% Viktor Witkovsky (witkovsky@savba.sk)
% Ver.: 17-Feb-2024 13:16:08

%% CHECK THE INPUTS AND OUTPUTS

if iscell(mu)
    mu = mu{1};
end

% Set the number of constraints
q = 9;
fun = zeros(q,1);

%% Constraints based on equation (8) in [Zuda (2023)] / Environment 1

fun(1) = mu(1)*mu(16) - ((1+beta(1))*(1-mu(12)*mu(14))) ...
                            + ((1+mu(10))*(1-mu(12)*mu(13))) ...
                            + (mu(21)*(1+mu(17)*mu(18))*(beta(3)-mu(11)));

fun(2) = mu(2)*mu(16) - ((1+beta(2))*(1-mu(12)*mu(15))) ...
                            + ((1+mu(10))*(1-mu(12)*mu(13))) ...
                            + (mu(21)*(1+mu(17)*mu(18))*(beta(4)-mu(11)));

fun(3) = mu(3)*mu(16) - ((1+beta(2))*(1-mu(12)*mu(15))) ...
                            + ((1+beta(1))*(1-mu(12)*mu(14))) ...
                            + (mu(21)*(1+mu(17)*mu(18))*(beta(4)-beta(3)));

%% Constraints based on equation (8) in [Zuda (2023)] / Environment 2

fun(4) = mu(4)*mu(16) - ((1+beta(1))*(1-mu(12)*mu(14))) ...
                            + ((1+mu(10))*(1-mu(12)*mu(13))) ...
                            + (mu(22)*(1+mu(17)*mu(19))*(beta(3)-mu(11)));

fun(5) = mu(5)*mu(16) - ((1+beta(2))*(1-mu(12)*mu(15))) ...
                            + ((1+mu(10))*(1-mu(12)*mu(13))) ...
                            + (mu(22)*(1+mu(17)*mu(19))*(beta(4)-mu(11)));

fun(6) = mu(6)*mu(16) - ((1+beta(2))*(1-mu(12)*mu(15))) ...
                            + ((1+beta(1))*(1-mu(12)*mu(14))) ...
                            + (mu(22)*(1+mu(17)*mu(19))*(beta(4)-beta(3)));

%% Constraints based on equation (8) in [Zuda (2023)] / Environment 3

fun(7) = mu(7)*mu(16) - ((1+beta(1))*(1-mu(12)*mu(14))) ...
                            + ((1+mu(10))*(1-mu(12)*mu(13))) ...
                            + (mu(23)*(1+mu(17)*mu(20))*(beta(3)-mu(11)));

fun(8) = mu(8)*mu(16) - ((1+beta(2))*(1-mu(12)*mu(15))) ...
                            + ((1+mu(10))*(1-mu(12)*mu(13))) ...
                            + (mu(23)*(1+mu(17)*mu(20))*(beta(4)-mu(11)));

fun(9) = mu(9)*mu(16) - ((1+beta(2))*(1-mu(12)*mu(15))) ...
                            + ((1+beta(1))*(1-mu(12)*mu(14))) ...
                            + (mu(23)*(1+mu(17)*mu(20))*(beta(4)-beta(3)));

end