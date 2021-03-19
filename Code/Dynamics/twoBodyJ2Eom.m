function dX = twoBodyJ2Eom(t,X,mu,J2,Rbody,wBody,fBody)
%% Local variables
r_i = X(1:3);
rMag = norm(r);

%% Higher Order Gravity
alpha = wBody*t;
T_e2i = [cos(alpha),-sin(alpha),0;...
         sin(alpha),cos(alpha),0;...
         0,0,1];
T_i2e = T_e2i';
r_e   = T_i2e*r_i;
lla = ecef2lla(r_e,fBody,Rbody);
lat = lla(1);

% J2 force potential
U_j2 = -mu*Rbody^2*J2*(3*sin(lat)^2-1)/2/rMag^3;

% partials wrt J2 force potential
dUdr   = 3*mu*Rbody^2*J2*(2*sin(lat)^2-1)/2/rMag^4;
dUdphi = -3*mu*Rbody^2*J2*sin(lat)*cos(lat)/rMag^3;

fJ2_geo = [dUdr;(1/rMag)*dUdphi;0];

fJ2_e = lla2ecef(fJ2_geo,fBody,Rbody);
fJ2_i = T_e2i*fJ2_e;
%% Derivatives
dX      = zeros(size(X));
dX(1:3) = X(4:6);
dX(4:6) = -mu*r/rMag+fJ2_i;