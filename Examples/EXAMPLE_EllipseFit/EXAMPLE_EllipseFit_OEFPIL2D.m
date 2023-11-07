%% EXAMPLE_EllipseFit_OEFPIL2D
%  Ellipse fit by OEFPIL2D / Artificial data

clear
close all

%% EXAMPLE (Ellipse fit by OEFPIL2D)

 x     = [ 0.9923    0.9820    0.8156    0.6528    0.5613 ...
           0.1245   -0.0835   -0.1754   -0.4786   -0.7933 ...
          -0.9044   -0.8952    0.0077    0.1180    0.2713]';

 y     = [ 0.6329    0.7410    0.7220    0.6971    0.6394 ...
           0.4553    0.2383    0.1550   -0.1082   -0.2388 ...
          -0.4654   -0.7143   -0.4161   -0.3032   -0.2883]';

 q     = 15;
 Ux    = 0.05^2*eye(q);
 Uy    = 0.05^2*eye(q);
 Uxy   = zeros(q);
 U     = [Ux Uxy; Uxy' Uy];
 fun   = @(mu,nu,beta) mu.^2 + beta(1)*nu.^2 + beta(2)*mu.*nu + ...
         beta(3)*mu + beta(4)*nu + beta(5);
 mu0   = x;
 nu0   = y;
 beta0 = [1; 0; 0; 0; -0.1];
 options.method = 'oefpil';
 result2D = OEFPIL2D(x,y,U,fun,mu0,nu0,beta0,options);