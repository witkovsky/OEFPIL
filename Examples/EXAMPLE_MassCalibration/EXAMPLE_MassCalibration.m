%% EXAMPLE / CALIBRATION of WEIGHT MASSES
%
% Viktor Witkovsky (witkovsky@savba.sk)
% Ver.: '14-Feb-2024 17:57:36'

clear
close all

%% CALIBRATION of WEIGHT MASSES
%
% Problem specified by Jaroslav ZUDA / Czech Metrology Institute (CMI) 
% E-MAIL: jzuda@cmi.cz
% Related EMPIR Project: 
% Improvement of the realisation of the mass scale 
% Short Name: RealMass, Project Number: 19RPT02
% 
% The fundamental challenge faced by any primary laboratory is ensuring the
% accurate calibration of weights. While the procedure for routine
% calibration of weights of lower accuracy classes is well-documented in
% OIML R111, procedures for calibration at the highest level are not as
% extensively covered. This issue is addressed by the RealMass project,
% which aims to develop such procedures.
% 
% Calibration at the highest level often requires the use of the split
% method for weight units. This procedure involves deriving individual
% weights from one or more reference weights using a measurement system
% with known weight values. Typically, these reference weights are
% prototypes made of a Pt and Ir alloy, although some institutes may use
% stainless steel weights. The weight is then calculated using the least
% squares method.
% 
% A significant challenge in accurate weight measurement is accounting for
% buoyancy forces. To address this, it is necessary to determine the
% density of the surrounding medium, which affects either the density or
% the volume of the weights. Additionally, the temperature of the
% environment must be considered since the volume of the weights changes
% with temperature.
% 
% There are various methods for determining the volume of weights, such as
% measuring the apparent weight in a liquid or comparing the apparent
% masses with a known volume body in the liquid. In recent years, there has
% been progress in measuring weights in a vacuum, particularly in light of
% the new definition of the unit of weight. However, measuring in a vacuum
% presents challenges due to surface phenomena on the weights. Therefore,
% it is often more practical to measure at low pressure, approximately
% corresponding to 50% atmospheric pressure, where surface changes are
% negligible.
% 
% Measurements in air and at low pressure can effectively determine both
% the volume and weight of weights. Comparisons indicate that this method
% yields results comparable to determining volume in a liquid.
% Consequently, there is interest in whether low-pressure measurement can
% be utilized for current weight and volume determination using the split
% method for weight units.

%% Load data

load EXAMPLE_MassCalibration.mat

%% Function of implicite restrictions on the model parameters

% function fun = funMassCalibration(mu,beta)
% % funMassCalibration - defines the necessary restrictions on the EIV model
% %  parameters for the calibration of weight masses. 
% % 
% %  funMassCalibration is an anonymous function that takes arguments mu and
% %  beta and returns fun, a q-dimensional column vector of constraints.
% %  Ideally, these constraints should equal zero. In this context, mu can be
% %  specified as either an m-dimensional column vector or as a cell array
% %  with one column vector {mu}, while beta is a p-dimensional vector. The
% %  output fun from funMassCalibration(mu, beta) represents the q
% %  constraints on the model parameters mu and beta. If the model
% %  restrictions are satisfied, fun = funMassCalibration(mu, beta) = 0.
% %   
% %  Variable related to the vector parameter mu:
% %
% %  mu(1)  = X1 = X1R_1  % comparison measurement: M1 - MR in evironment 1
% %  mu(2)  = X2 = X2R_1  % comparison measurement: M2 - MR in evironment 1
% %  mu(3)  = X3 = X21_1  % comparison measurement: M2 - M1 in evironment 1
% %  mu(4)  = X4 = X1R_2  % comparison measurement: M1 - MR in evironment 2
% %  mu(5)  = X5 = X2R_2  % comparison measurement: M2 - MR in evironment 2
% %  mu(6)  = X6 = X21_2  % comparison measurement: M2 - M1 in evironment 2
% %  mu(7)  = X7 = X1R_3  % comparison measurement: M1 - MR in evironment 3
% %  mu(8)  = X8 = X2R_3  % comparison measurement: M2 - MR in evironment 3
% %  mu(9)  = X9 = X21_3  % comparison measurement: M2 - M1 in evironment 3
% %  mu(10) = dMR         % difference from nominal of the reference weight
% %  mu(11) = VR          % volume of the reference weight
% %  mu(12) = G           % gravitational ratio
% %  mu(13) = h1 = hR     % center of gravity of the reference weight
% %  mu(14) = h2 = h1     % center of gravity of the measured weight 1
% %  mu(15) = h3 = h2     % center of gravity of the measured weight 2
% %  mu(16) = K           % adjustment constant
% %  mu(17) = α           % coefficient of thermal expansion
% %  mu(18) = dT1 = dT_1  % temperature differnce T_1 - 20 in evironment 1
% %  mu(19) = dT2 = dT_2  % temperature differnce T_2 - 20 in evironment 2
% %  mu(20) = dT3 = dT_3  % temperature differnce T_3 - 20 in evironment 3
% %  mu(21) = ρ1  = ρ_1   % density of air in evironment 1
% %  mu(22) = ρ2  = ρ_2   % density of air in evironment 2
% %  mu(23) = ρ3  = ρ_3   % density of air in evironment 3
% %
% %  Variable related to the vector parameter beta:
% %
% %  beta(1) = dM1        % difference from nominal of the measured weight 1
% %  beta(2) = dM2        % difference from nominal of the measured weight 2
% %  beta(3) = V1         % volume of the measured weight 1
% %  beta(4) = V2         % volume of the measured weight 2
% %
% % EXAMPLE
% %  data  = [-9.7550e-08 -9.0350e-08  6.3200e-09 -4.7849e-07 9.4810e-08 ...
% %            3.8381e-07 -3.2604e-07 -9.3040e-08  2.3331e-07 6.4100e-07 ...
% %            1.2544e-04  3.0000e-07  3.6500e-02  3.6500e-02 3.0000e-02 ...
% %            9.9986e-01  4.8000e-05  8.0000e-01  4.0000e-01 9.0000e-01 ...
% %            1.1552e+00  5.8770e-01  8.1380e-01]';
% %  ux    = [ 5.0000e-10  3.0000e-10  2.0000e-10  3.0000e-10 3.0000e-10 ...  
% %            3.0000e-10  2.0000e-10  3.0000e-10  2.0000e-10 4.0000e-08 ...
% %            1.0000e-09  3.0000e-08  5.0000e-04  5.0000e-04 5.0000e-04 ... 
% %            5.0000e-06  2.0000e-06  1.0000e-01  1.0000e-01 1.0000e-01 ...
% %            5.0000e-04  5.0000e-04  5.0000e-04]';
% %  U     = diag(ux.^2);
% %  fun   = @(mu,beta) funMassCalibration(mu,beta)
% %  mu0   = data;
% %  beta0 = [mu0(10) mu0(11) mu0(10) mu0(11)]';
% %  clear options
% %  options.method        = 'oefpilrs2';
% %  options.criterion     = 'function';
% %  options.numDiffMethod = 'LynnesMoller';
% %  options.q             = 9;
% %  options.tol           = 1e-15;
% %  result = OEFPIL(data,U,fun,mu0,beta0,options);
% 
% % Viktor Witkovsky (witkovsky@savba.sk)
% % Ver.: 17-Feb-2024 13:16:08
% 
% %% CHECK THE INPUTS AND OUTPUTS
% 
% if iscell(mu)
%     mu = mu{1};
% end
% 
% % Set the number of constraints
% q = 9;
% fun = zeros(q,1);
% 
% %% Constraints based on equation (8) in [Zuda (2023)] / Environment 1
% 
% fun(1) = mu(1)*mu(16) - ((1+beta(1))*(1-mu(12)*mu(14))) ...
%                             + ((1+mu(10))*(1-mu(12)*mu(13))) ...
%                             + (mu(21)*(1+mu(17)*mu(18))*(beta(3)-mu(11)));
% 
% fun(2) = mu(2)*mu(16) - ((1+beta(2))*(1-mu(12)*mu(15))) ...
%                             + ((1+mu(10))*(1-mu(12)*mu(13))) ...
%                             + (mu(21)*(1+mu(17)*mu(18))*(beta(4)-mu(11)));
% 
% fun(3) = mu(3)*mu(16) - ((1+beta(2))*(1-mu(12)*mu(15))) ...
%                             + ((1+beta(1))*(1-mu(12)*mu(14))) ...
%                             + (mu(21)*(1+mu(17)*mu(18))*(beta(4)-beta(3)));
% 
% %% Constraints based on equation (8) in [Zuda (2023)] / Environment 2
% 
% fun(4) = mu(4)*mu(16) - ((1+beta(1))*(1-mu(12)*mu(14))) ...
%                             + ((1+mu(10))*(1-mu(12)*mu(13))) ...
%                             + (mu(22)*(1+mu(17)*mu(19))*(beta(3)-mu(11)));
% 
% fun(5) = mu(5)*mu(16) - ((1+beta(2))*(1-mu(12)*mu(15))) ...
%                             + ((1+mu(10))*(1-mu(12)*mu(13))) ...
%                             + (mu(22)*(1+mu(17)*mu(19))*(beta(4)-mu(11)));
% 
% fun(6) = mu(6)*mu(16) - ((1+beta(2))*(1-mu(12)*mu(15))) ...
%                             + ((1+beta(1))*(1-mu(12)*mu(14))) ...
%                             + (mu(22)*(1+mu(17)*mu(19))*(beta(4)-beta(3)));
% 
% %% Constraints based on equation (8) in [Zuda (2023)] / Environment 3
% 
% fun(7) = mu(7)*mu(16) - ((1+beta(1))*(1-mu(12)*mu(14))) ...
%                             + ((1+mu(10))*(1-mu(12)*mu(13))) ...
%                             + (mu(23)*(1+mu(17)*mu(20))*(beta(3)-mu(11)));
% 
% fun(8) = mu(8)*mu(16) - ((1+beta(2))*(1-mu(12)*mu(15))) ...
%                             + ((1+mu(10))*(1-mu(12)*mu(13))) ...
%                             + (mu(23)*(1+mu(17)*mu(20))*(beta(4)-mu(11)));
% 
% fun(9) = mu(9)*mu(16) - ((1+beta(2))*(1-mu(12)*mu(15))) ...
%                             + ((1+beta(1))*(1-mu(12)*mu(14))) ...
%                             + (mu(23)*(1+mu(17)*mu(20))*(beta(4)-beta(3)));
% 
% end

%% EXAMPLE

%U     = U_diag;
%U     = U_full;
U     = U_full2;
fun   = @(mu,beta) funMassCalibration(mu,beta);
mu0   = data;
% beta0 = [mu0(10) mu0(11) mu0(10) mu0(11)]';
beta0 = [0 0 0 0]';

clear options
options.method        = 'oefpilrs2';
options.criterion     = 'function';
options.numDiffMethod = 'LynnesMoller';
options.q             = 9;
options.tol           = 1e-15;
options.isEstimatedVariance = false;

%% OEFPIL fit

result = OEFPIL(data,U,fun,mu0,beta0,options);

%% RESULTS

mufit = result.mu;

% Estimated (fitted) differences from nominal weight 1kg
dMRfit_51699 = mufit(10);       % Reference weight: 51699
dM1fit_15N   = result.beta(1);  % Measured weight: 15N
dM2fit_15963 = result.beta(2);  % Measured weight: 15936

dmass = [ dMRfit_51699;  dM1fit_15N; dM2fit_15963];
massNominal = 1;
mass = massNominal + dmass;

% Standard uncertainty of the weight
u_MRfit_51699 = result.umu(10);  % Uncertainty of the reference weight: 51699
u_M1fit_15N   = result.ubeta(1); % Uncertainty of the measured weight: 15N
u_M2fit_15963 = result.ubeta(2); % Uncertainty of the measured weight: 15936

u_mass = [ u_MRfit_51699;  u_M1fit_15N; u_M2fit_15963];

% Estimated (fitted) volumes of the weights
VRfit_51699 = mufit(11);        % Volume of the reference weight: 51699
V1fit_15N   = result.beta(3);   % Volume of the measured weight: 15N
V2fit_15963 = result.beta(4);   % Volume of the measured weight: 15936

volume = [ VRfit_51699;  V1fit_15N; V2fit_15963];

% Standard uncertainty of the volumes
u_VRfit_51699 = result.umu(11);  % Uncertainty of the volume of the reference weight: 51699
u_V1fit_15N   = result.ubeta(3); % Uncertainty of the volume of the measured weight: 15N
u_V2fit_15963 = result.ubeta(4); % Uncertainty of the volume of the measured weight: 15936

u_volume = [ u_VRfit_51699;  u_V1fit_15N; u_V2fit_15963];

%% Best Estimates by OEFPIL
BestEstimates   = [mass;volume];

disp(table(BestEstimates))

%    BestEstimates
% ____________________
%
%          1.000000641
%    0.999999887200819
%     1.00000075427854
%           0.00012544
% 0.000124891615012868
% 0.000125642839434949


%% Uncertainty (covariance) matrix of the OEFPIL Best Estimates 
%  Estimated parameters of interest (M_R, M_1, M_2, V_R, V_1, V_2):

id_BestEstimates = [10 24 25 11 26 27];
U_BestEstimates = result.Umubeta(id_BestEstimates,id_BestEstimates);

disp(U_BestEstimates)

% U_BestEstimates =
% 
%    1.0e-14 *
% 
%    0.159999999999415   0.159999999999991   0.159999999687991   0.000116400000000   0.000116400000000   0.000116400000000
%    0.159999999999991   0.160046483572160   0.160017224314753   0.000116400000000   0.000156928900609   0.000137367876887
%    0.159999999687991   0.160017224314753   0.160038988503685   0.000116399999773   0.000137232205722   0.000149667266298
%    0.000116400000000   0.000116400000000   0.000116399999773   0.000100000000000   0.000100000000000   0.000100000000000
%    0.000116400000000   0.000156928900609   0.000137232205722   0.000100000000000   0.000149348606540   0.000126375525396
%    0.000116400000000   0.000137367876887   0.000149667266298   0.000100000000000   0.000126375525396   0.000139639048366


%% Standard uncertainties of the OEFPIL Best Estimates 
%  Estimated parameters of interest (M_R, M_1, M_2, V_R, V_1, V_2):

u_BestEstimates = sqrt(diag(U_BestEstimates));

disp(table(u_BestEstimates))

%   u_BestEstimates
% ____________________
%
% 3.99999999999269e-08
% 4.00058100245652e-08
% 4.00048732661016e-08
%                1e-09
% 1.22208267535257e-09
% 1.18168967316397e-09




%% Correlation matrix of the OEFPIL Best Estimates 
%  Estimated parameters of interest (M_R, M_1, M_2, V_R, V_1, V_2):

Corr_BestEstimates = U_BestEstimates ./ (u_BestEstimates * u_BestEstimates');

disp(Corr_BestEstimates)

% Corr_BestEstimates =
% 
%    1.000000000000000   0.999854770482318   0.999878181240613   0.029100000000021   0.023811809615624   0.024625754680711
%    0.999854770482318   1.000000000000000   0.999840594581596   0.029095773820982   0.032098096361323   0.029057528836698
%    0.999878181240613   0.999840594581596   1.000000000000000   0.029096455074049   0.028070009357776   0.031659969132984
%    0.029100000000021   0.029095773820982   0.029096455074049   1.000000000000000   0.818275244521817   0.846245865314624
%    0.023811809615624   0.032098096361323   0.028070009357776   0.818275244521817   1.000000000000000   0.875102544210402
%    0.024625754680711   0.029057528836698   0.031659969132984   0.846245865314624   0.875102544210402   1.000000000000000

%% TABLE  OEFPIL Best Estimates of Mass
%  Estimated parameters of interest (M_R, M_1, M_2):

alpha = 0.05;
coverageFactor = norminv(1-alpha/2);

TABLE_mass= table;
TABLE_mass.Properties.Description = 'Estimated Mass of Weights by OEFPIL';
TABLE_mass.ESTIMATE = mass;
TABLE_mass.STD      = u_mass;
TABLE_mass.FACTOR   = coverageFactor*ones(size(mass));
TABLE_mass.LOWER    = mass - coverageFactor*u_mass;
TABLE_mass.UPPER    = mass + coverageFactor*u_mass;
TABLE_mass.PVAL     = 2*normcdf(-abs((1-mass)./u_mass));
TABLE_mass.Properties.RowNames = {'M_R (51699)' 'M_1 (15N)' 'M_2 (15963)'};
TABLE_mass.Properties.VariableNames = {'ESTIMATE [kg]'  'STD [kg]'  'FACTOR'  'LOWER BOUND [kg]'  'UPPER BOUND [kg]'  'PVAL [H0: Mass = 1 kg]'};

disp(TABLE_mass)

%                  ESTIMATE [kg]            STD [kg]               FACTOR         LOWER BOUND [kg]     UPPER BOUND [kg]     PVAL [H0: Mass = 1 kg]
%                _________________    ____________________    ________________    _________________    _________________    ______________________
%
% M_R (51699)          1.000000641    3.99999999999269e-08    1.95996398454005     1.00000056260144     1.00000071939856      8.5490885926873e-58
% M_1 (15N)      0.999999887200819    4.00058100245652e-08    1.95996398454005    0.999999808790873    0.999999965610766      0.00480880457852863
% M_2 (15963)     1.00000075427854    4.00048732661016e-08    1.95996398454005     1.00000067587043     1.00000083268665     2.69030217404407e-79

%% TABLE  OEFPIL Best Estimates of Volume
%  Estimated parameters of interest (V_R, V_1, V_2):

TABLE_volume= table;
TABLE_volume.Properties.Description = 'Estimated Volume of Weights by OEFPIL';
TABLE_volume.ESTIMATE = volume;
TABLE_volume.STD      = u_volume;
TABLE_volume.FACTOR   = coverageFactor*ones(size(volume));
TABLE_volume.LOWER    = volume - coverageFactor*u_volume;
TABLE_volume.UPPER    = volume + coverageFactor*u_volume;
TABLE_volume.PVAL     = 2*normcdf(-abs((0.000125-volume)./u_volume));
TABLE_volume.Properties.RowNames = {'V_R (51699)' 'V_1 (15N)' 'V_2 (15963)'};
TABLE_volume.Properties.VariableNames = {'ESTIMATE [m^3]'  'STD [m^3]'  'FACTOR'  'LOWER BOUND [m^3]'  'UPPER BOUND [m^3]'  'PVAL [H0: Volume = 125 cm^3]'};

disp(TABLE_volume)

%                   ESTIMATE [m^3]            STD [m^3]               FACTOR          LOWER BOUND [m^3]       UPPER BOUND [m^3]      PVAL [H0: Volume = 125 cm^3]
%                ____________________    ____________________    ________________    ____________________    ____________________    ____________________________
%
% V_R (51699)              0.00012544                   1e-09    1.95996398454005    0.000125438040036015    0.000125441959963985                 0
% V_1 (15N)      0.000124891615012868    1.22208267535257e-09    1.95996398454005    0.000124889219774838    0.000124894010250898                 0
% V_2 (15963)    0.000125642839434949    1.18168967316397e-09    1.95996398454005    0.000125640523365748    0.000125645155504149                 0

%% TABLE / Estimated parameters / OEFPIL

p = length(BestEstimates);
TABLE_OEFPIL = table;
TABLE_OEFPIL.Properties.Description = 'Estimated Mass and Volume of the Weights by OEFPIL';
TABLE_OEFPIL.ESTIMATE = BestEstimates;
TABLE_OEFPIL.STD      = u_BestEstimates;
TABLE_OEFPIL.FACTOR   = coverageFactor*ones(size(BestEstimates));
TABLE_OEFPIL.LOWER    = BestEstimates - coverageFactor*u_BestEstimates;
TABLE_OEFPIL.UPPER    = BestEstimates + coverageFactor*u_BestEstimates;
TABLE_OEFPIL.Properties.RowNames = {'M_R (51699)' 'M_1 (15N)' 'M_2 (15963)' 'V_R (51699)' 'V_1 (15N)' 'V_2 (15963)'};
TABLE_OEFPIL.Properties.VariableNames = {'ESTIMATE [kg | m^3]'  'STD [kg | m^3]'  'FACTOR'  'LOWER BOUND [kg | m^3]'  'UPPER BOUND [kg | m^3]' };

disp(TABLE_OEFPIL)

%                ESTIMATE [kg | m^3]        STD [kg | m^3]            FACTOR         LOWER BOUND [kg | m^3]    UPPER BOUND [kg | m^3]
%                ____________________    ____________________    ________________    ______________________    ______________________
%
% M_R (51699)             1.000000641    3.99999999999269e-08    1.95996398454005         1.00000056260144          1.00000071939856
% M_1 (15N)         0.999999887200819    4.00058100245652e-08    1.95996398454005        0.999999808790873         0.999999965610766
% M_2 (15963)        1.00000075427854    4.00048732661016e-08    1.95996398454005         1.00000067587043          1.00000083268665
% V_R (51699)              0.00012544                   1e-09    1.95996398454005     0.000125438040036015      0.000125441959963985
% V_1 (15N)      0.000124891615012868    1.22208267535257e-09    1.95996398454005     0.000124889219774838      0.000124894010250898
% V_2 (15963)    0.000125642839434949    1.18168967316397e-09    1.95996398454005     0.000125640523365748      0.000125645155504149


%% Save the results

save RESULTS_MassCalibration.mat