function dX = twoBodyEom(t,X,vProc)
%% Local variables
r = X(1:3);
rMag = norm(r);
mu = X(7);

%% Derivatives
dX      = zeros(size(X));
dX(1:3) = X(4:6);
dX(4:6) = -mu*r/rMag^3+vProc;
dX(7)   = 0;