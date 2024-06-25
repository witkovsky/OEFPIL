%% EXAMPLE_RobotArmSetting_DataEllipse
%
% This example is inspired by the paper [1], submitted for publication in
% Measurement Science Review by Krepelka, Schovanek, Tucek, Hrabovsky, and
% Jane in 2023 (accepted for publication in 2024).
%
% Here we use the artificial 3D ellipse data presented in [1] (see Table 2).

% Viktor Witkovsky (witkovsky@savba.sk)
% Ver.: 07-Jan-2024 12:29:56

clear 
close all

%% Load Ellipse data

load DataEllipse.mat

%% Ellipse Data - Table 2 from the paper Krepelka etal MSR 2024

% Measured data
XYZdata = [ ...
  1358.874                  -250.395                   660.817
  1359.712                  -212.453                   646.804
  1337.276                   -172.65                   621.252
  1297.452                  -131.235                   584.957
  1237.693                    -84.54                   542.004
  1160.935                   -35.854                   490.329
  1064.542                    11.942                   432.137
   949.067                    61.656                   366.842
   822.277                   110.383                   296.436
   681.333                   155.583                   221.627
   528.406                   204.025                   141.752
   369.478                   245.877                    67.816
   204.851                   284.963                   -13.342
     39.42                   319.494                   -92.695
  -129.532                   353.129                  -171.121
  -293.352                   378.851                  -245.085
  -452.946                    397.15                  -314.136
  -606.428                   412.267                  -378.611
  -748.978                   421.462                  -435.984
  -880.648                   424.912                  -484.202
  -993.653                   421.318                  -527.011
 -1095.357                   410.078                  -561.702
 -1177.098                   396.412                  -584.088
 -1241.456                   372.772                   -598.56
 -1285.651                   349.035                  -601.591
 -1305.739                   316.665                  -597.842
 -1309.606                   279.743                  -579.359
 -1292.067                   241.963                  -555.232
 -1250.883                   196.316                  -520.326
 -1190.357                   152.203                  -476.362
 -1110.784                   104.508                   -424.75
 -1014.365                    56.242                  -364.004
  -901.891                     6.167                  -301.974
  -773.813                   -42.472                  -229.988
  -631.703                   -88.264                  -153.169
  -482.454                  -135.015                   -78.914
  -324.489                  -178.371                     0.812
  -158.032                  -217.965                    79.767
     8.905                  -253.564                   160.522
   177.637                  -283.075                   236.949
   342.537                  -309.395                   311.713
   503.871                  -330.763                   379.369
   654.951                   -345.77                   442.161
   798.591                  -353.444                    502.72
   927.696                  -356.569                   552.527
  1042.979                  -352.494                   592.797
  1144.224                  -342.742                   627.773
  1225.422                  -327.735                    650.25
  1291.302                  -306.727                   664.902
  1333.763                  -280.449                   667.358];

% Target position
Target = [ ...
  1419.539                        20                    528.03
  1410.657                    57.029                   517.849
  1379.687                    93.474                   499.738
  1327.115                    128.76                   473.983
  1253.772                    162.33                   440.989
  1160.814                   193.657                   401.276
  1049.706                   222.244                   355.472
   922.202                   247.642                   304.299
   780.312                    269.45                   248.563
   626.274                   287.324                   189.144
   462.517                   300.982                   126.978
   291.623                   310.209                    63.047
   116.288                   314.859                    -1.643
   -60.723                   314.859                    -66.07
  -236.619                   310.209                  -129.218
  -408.626                   300.982                  -190.092
  -574.031                   287.324                  -247.731
  -730.225                    269.45                  -301.228
  -874.745                   247.642                  -349.737
 -1005.313                   222.244                  -392.493
 -1119.868                   193.657                  -428.824
 -1216.605                    162.33                  -458.155
 -1293.997                    128.76                  -480.024
 -1350.825                    93.474                  -494.086
 -1386.191                    57.029                   -500.12
 -1399.539                        20                   -498.03
 -1390.657                   -17.029                  -487.849
 -1359.687                   -53.474                  -469.738
 -1307.115                    -88.76                  -443.983
 -1233.772                   -122.33                  -410.989
 -1140.814                  -153.657                  -371.276
 -1029.706                  -182.244                  -325.472
  -902.202                  -207.642                  -274.299
  -760.312                   -229.45                  -218.563
  -606.274                  -247.324                  -159.144
  -442.517                  -260.982                   -96.978
  -271.623                  -270.209                   -33.047
   -96.288                  -274.859                    31.643
    80.723                  -274.859                     96.07
   256.619                  -270.209                   159.218
   428.626                  -260.982                   220.092
   594.031                  -247.324                   277.731
   750.225                   -229.45                   331.228
   894.745                  -207.642                   379.737
  1025.313                  -182.244                   422.493
  1139.868                  -153.657                   458.824
  1236.605                   -122.33                   488.155
  1313.997                    -88.76                   510.024
  1370.825                   -53.474                   524.086
  1406.191                   -17.029                    530.12];

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
% **corrected** points and in Table 2 are presented the following results:

% Corrected = [ ...
%       1422.341                    19.381                   526.316
%       1412.157                    57.722                    518.38
%       1378.091                    93.203                   500.461
%       1325.937                   127.051                   473.024
%        1252.66                    161.52                   440.793
%        1161.75                   194.543                   401.078
%       1051.605                   222.363                    355.73
%        921.982                   247.887                   304.584
%        781.279                   270.147                    248.78
%        627.478                   285.866                   188.756
%        460.868                   302.168                   124.812
%        290.699                    309.76                    65.983
%        115.127                   314.187                    -0.132
%        -59.889                   313.806                   -65.139
%       -237.982                     311.5                  -129.204
%        -408.71                   302.406                  -190.348
%       -573.031                   286.686                  -248.007
%       -730.253                   268.934                   -301.94
%       -874.808                    247.48                  -350.313
%      -1006.617                   222.195                  -391.069
%      -1118.178                   194.441                   -428.48
%      -1216.186                   161.298                  -459.615
%      -1293.285                   129.485                   -479.94
%      -1350.397                    91.719                  -494.846
%      -1387.083                    57.676                  -499.439
%      -1397.755                    20.915                  -499.842
%      -1390.359                   -17.881                  -487.194
%      -1361.649                   -52.546                  -470.131
%      -1307.219                   -89.875                  -444.789
%       -1233.78                  -121.815                  -411.194
%      -1140.421                  -153.131                  -371.521
%      -1029.907                  -181.745                  -323.728
%       -903.364                  -207.966                  -275.717
%       -761.306                  -230.016                  -218.402
%       -606.055                  -246.264                  -156.543
%       -443.912                  -261.161                   -97.704
%       -273.776                  -271.245                   -33.352
%        -96.504                  -275.362                    30.423
%         80.358                  -275.791                    96.575
%        257.004                  -269.335                   159.323
%        428.994                  -260.633                   221.112
%        595.659                  -247.181                   276.832
%        750.393                  -229.735                   329.235
%        895.884                  -206.979                   380.937
%       1025.279                  -182.291                   423.434
%       1138.687                  -153.116                   458.307
%       1236.662                  -121.556                   489.513
%       1312.901                   -88.621                   510.164
%       1372.171                   -52.982                   524.762
%       1406.483                   -16.901                   529.269];

% For comparison, see the transformed OEFPIL data (rounded to three decimal
% places) given in TransformedData 
% TransformedData = [ ...
%       1422.344                    19.374                   526.306
%       1412.161                    57.715                   518.372
%       1378.093                    93.196                   500.453
%        1325.94                   127.045                   473.017
%       1252.662                   161.514                   440.788
%       1161.751                   194.537                   401.075
%       1051.605                   222.358                   355.729
%        921.982                   247.882                   304.584
%        781.278                   270.143                   248.782
%        627.476                   285.861                    188.76
%        460.866                   302.164                   124.817
%        290.696                   309.756                     65.99
%        115.123                   314.184                    -0.124
%        -59.893                   313.803                    -65.13
%       -237.986                   311.498                  -129.194
%       -408.714                   302.404                  -190.337
%       -573.036                   286.684                  -247.995
%       -730.257                   268.933                  -301.927
%       -874.814                   247.479                  -350.299
%      -1006.622                   222.194                  -391.056
%      -1118.183                    194.44                  -428.467
%      -1216.191                   161.297                  -459.602
%       -1293.29                   129.484                  -479.928
%      -1350.402                    91.718                  -494.834
%      -1387.087                    57.675                  -499.428
%      -1397.758                    20.914                  -499.832
%      -1390.362                   -17.882                  -487.184
%      -1361.652                   -52.548                  -470.123
%      -1307.222                   -89.877                  -444.783
%      -1233.782                  -121.818                  -411.189
%      -1140.423                  -153.134                  -371.517
%      -1029.908                  -181.749                  -323.727
%       -903.365                  -207.969                  -275.717
%       -761.306                   -230.02                  -218.404
%       -606.055                  -246.268                  -156.546
%        -443.91                  -261.166                   -97.708
%       -273.774                  -271.251                   -33.358
%        -96.502                  -275.368                    30.416
%         80.361                  -275.797                    96.566
%        257.007                  -269.342                   159.313
%        428.998                   -260.64                   221.101
%        595.662                  -247.188                    276.82
%        750.397                  -229.742                   329.223
%        895.889                  -206.987                   380.924
%       1025.284                  -182.298                   423.421
%       1138.691                  -153.124                   458.294
%       1236.666                  -121.563                     489.5
%       1312.905                   -88.628                   510.151
%       1372.175                    -52.99                   524.751
%       1406.487                   -16.908                   529.258];

% Here is the difference between the OEFPIL and original solution in [1]:

Differences = TransformedData - Corrected;
Tab6 = table(Differences);
disp(Tab6)

