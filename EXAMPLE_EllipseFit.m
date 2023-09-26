
% EXAMPLE 1: (Fit the ellipse for generated measurements x and y)
% D:\VW\VW Odborne\VW Publications\VW-2013\VW-2013-h (Koening)\Matlab
%
alpha0 = 0; beta0 = 0;        % true ellipse center [0,0]
alpha1 = 1; beta1 = 0.75;     % true amplitudes
phi0 = pi/3;                  % phase shift
X = @(t) alpha0 + alpha1 * cos(t);
Y = @(t) beta0 + beta1 * sin(t + phi0);
sigma = 0.05;                 % true error STD
N = 15;                       % No. of observations (x,y)
Ncycles = 0.8;                % No. of whole ellipse cycles
phi = Ncycles*(2*pi)*sort(rand(N,1)); % true phases
x = X(phi) + sigma*randn(size(X(phi)));
y = Y(phi) + sigma*randn(size(Y(phi)));
result2 = ellipseFit4HC(x,y);
disp(result.TABLE_ellipsePars)

Xf = result2.X_fitfun;
Yf = result2.Y_fitfun;
t = linspace(0,2*pi,101);
xx = Xf(t);
yy = Yf(t);
Fitted ellipse
plot(xx,yy,'-')

%% EXAMPLE (Ellipse fit by OEFPIL)
 x     = [ 1.0310    1.0150    0.9677    0.8924    0.8025 ...
          -0.3599   -0.7536   -0.7702   -1.0107   -0.9757 ...
          -0.8163   -0.6509   -0.3045   -0.0935    0.3155]';

 y     = [ 0.7332    0.6956    0.7028    0.7678    0.7955 ...
           0.1984   -0.3282   -0.3222   -0.7162   -0.7561 ...
          -0.6521   -0.7264   -0.6397   -0.4567   -0.2040]';
 m     = 15;
 Ux    = 0.05^2*eye(m);
 Uy    = 0.05^2*eye(m);
 Uxy   = zeros(m);
 U     = [Ux Uxy; Uxy' Uy];
 fun   = @(mu,nu,beta) mu.^2 + beta(1)*nu.^2 + beta(2)*mu.*nu + ...
         beta(3)*mu + beta(4)*nu + beta(5);
 mu0   = x;
 nu0   = y;
 beta0 = [1; 0; 0; 0; -0.1];
 options.method = 'oefpil';
 result = OEFPIL2D(x,y,U,fun,mu0,nu0,beta0,options);