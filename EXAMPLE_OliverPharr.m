%% EXAMPLE: Oliver-Phar function fit / fun(mu,nu,[0.75;-0.25;1.75]))

 x     = [ 0.2505    2.6846    5.1221    7.5628   10.0018 ...
           0.2565    2.6858    5.1255    7.5623    9.9952 ...
           0.2489    2.6830    5.1271    7.5603   10.0003]';

 y     = [ 0.2398    4.9412   14.2090   27.3720   44.0513 ...
           0.2110    4.9517   14.2306   27.3937   44.0172 ...
           0.2303    4.9406   14.2690   27.3982   44.0611]';

 m     = 5;
 Ux    = 0.05^2*eye(3*m);
 Uyblk = 0.1^2*eye(m)+ 0.05^2*ones(m);
 Uy    = blkdiag(Uyblk,Uyblk,Uyblk);
 Uxy   = zeros(3*m);
 U     = [Ux Uxy; Uxy' Uy];

 fun   = @(mu,nu,beta) beta(1).*(mu-beta(2)).^beta(3) - nu;
 
 mu0   = x;
 nu0   = y;
 beta0 = [1 0 2]';
 
 options.method = 'oefpil';
 options.funDiff_mu = @(mu,nu,beta) beta(1)*beta(3)*(mu - beta(2)).^(beta(3)-1);
 options.funDiff_nu = @(mu,nu,beta) - ones(size(nu));
 options.funDiff_beta = @(mu,nu,beta) [(mu-beta(2)).^beta(3), ...
                      -beta(1)*beta(3)*(mu-beta(2)).^(beta(3)-1), ...
                       beta(1).*(mu-beta(2)).^beta(3).*log(mu-beta(2))];
 
 result = OEFPIL2D(x,y,U,fun,mu0,nu0,beta0,options);