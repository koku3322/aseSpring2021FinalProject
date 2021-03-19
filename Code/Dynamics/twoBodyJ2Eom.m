function dX = twoBodyJ2Eom(t,X,mu,J2,Rbody,wBody)
%% Local variables
r_i = X(1:3);
rMag = norm(r_i);

%% Higher Order Gravity
alpha = wBody*t;
T_e2i = [cos(alpha),-sin(alpha),0;...
         sin(alpha),cos(alpha),0;...
         0,0,1];
T_i2e = T_e2i';
r_e   = T_i2e*r_i;
phi = acos(r_e(3)/rMag);%latitude
lam = atan2(r_e(2),r_e(1));%longitude

% partials wrt J2 force potential
dUdr   = 3*mu*Rbody^2*J2*(2*sin(phi)^2-1)/2/rMag^4;
dUdphi = -3*mu*Rbody^2*J2*sin(phi)*cos(phi)/rMag^3;

fJ2_geo = [dUdr;(1/rMag)*dUdphi;0];

T_g2e = [cos(phi)*cos(lam),-sin(phi)*cos(lam),-sin(lam);...
         cos(phi)*sin(lam),-sin(phi)*sin(lam),cos(lam);...
         sin(lam),cos(lam),0];
     
fJ2_e = T_g2e*fJ2_geo;
fJ2_i = T_e2i*fJ2_e;
%% Derivatives
dX      = zeros(size(X));
dX(1:3) = X(4:6);
dX(4:6) = -mu*r_i/rMag^3+fJ2_i;