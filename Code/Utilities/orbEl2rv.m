function [r,v] = orbEl2rv(a,ecc,inc,w,ra,f,mu)

% lcoal variables
cw = cos(w);cW = cos(ra);ci = cos(inc);
sw = sin(w);sW = sin(ra);si = sin(inc);

%% Compute position and velocity in the orbital frame
r = a*(1-ecc^2)/(1+ecc*cos(f));

rOrbit = r*[cos(f);sin(f);0];
vOrbit = sqrt(mu/a/(1-ecc^2))*[-sin(f);ecc+cos(f);0];

%% Convert to cartesian frame
rx = rOrbit(1)*( cw*cW - sw*ci*sW ) - rOrbit(2)*( sw*cW + cw*ci*sW );
ry = rOrbit(1)*( cw*sW + sw*ci*cW ) + rOrbit(2)*( cw*ci*cW - sw*sW );
rz = rOrbit(1)*( sw*si ) + rOrbit(2)*( cw*si );

r = [rx;ry;rz];

vx = vOrbit(1)*( cw*cW - sw*ci*sW ) - vOrbit(2)*( sw*cW + cw*ci*sW );
vy = vOrbit(1)*( cw*sW + sw*ci*cW ) + vOrbit(2)*( cw*ci*cW - sw*sW );
vz = vOrbit(1)*( sw*si ) + vOrbit(2)*( cw*si );

v = [vx;vy;vz];