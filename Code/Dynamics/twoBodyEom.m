function dX = twoBodyEom(t,X,mu,vProc)
%% Local variables
r = X(1:3);
rMag = norm(r);

%% Derivatives
dX      = zeros(size(X));
dX(1:3) = X(4:6);
dX(4:6) = -mu*r/rMag^3+vProc;