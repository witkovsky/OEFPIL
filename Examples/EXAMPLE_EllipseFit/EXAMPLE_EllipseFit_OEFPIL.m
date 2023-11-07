%% EXAMPLE_EllipseFit_OEFPIL
%  Ellipse fit by OEFPIL / Artificial data

clear
close all
 
 %% EXAMPLE (Ellipse fit by OEFPIL)

 x     = [ 0.9923    0.9820    0.8156    0.6528    0.5613 ...
           0.1245   -0.0835   -0.1754   -0.4786   -0.7933 ...
          -0.9044   -0.8952    0.0077    0.1180    0.2713]';

 y     = [ 0.6329    0.7410    0.7220    0.6971    0.6394 ...
           0.4553    0.2383    0.1550   -0.1082   -0.2388 ...
          -0.4654   -0.7143   -0.4161   -0.3032   -0.2883]';
 
 q     = length(x);
 Ux    = 0.05^2*eye(q);
 Uy    = 0.05^2*eye(q);
 fun   = @(mu,beta) mu{1}.^2 + beta(1)*mu{2}.^2 + beta(2)*mu{1}.*mu{2} ...
         + beta(3)*mu{1} + beta(4)*mu{2} + beta(5);
 beta0 = [1; 0; 0; 0; -0.1];
 options.method = 'oefpil';

 result = OEFPIL({x,y},{Ux,Uy},fun,{x,y},beta0,options);
 disp(result)