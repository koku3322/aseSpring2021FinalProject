function [a,ecc,inc,w,ra,f] = rvt2OrbEl(r0,v0,mu)
% This function calculates orbital elements from position, velocity and
% time

%% Define local variables
zHat = [0;0;1];
rMag = norm(r0);

%% Eccentricity
eVec = (1/mu)*cross(v0,cross(r0,v0)) - r0/rMag;

ecc = norm(eVec);

%% Semi-major axis
h = cross(r0,v0);
hMag = norm(h);

a = hMag^2/(mu*(1-ecc^2));

%% Inclination

inc = acos(h(3)/hMag); % radians

%% Right Ascension
n = cross(zHat,h);% node vector

if n(2) >=0
    ra = acos(n(1)/norm(n));
else
    ra = 2*pi-acos(n(1)/norm(n));
end
%% Arguement of Periapsis
if eVec(3) >= 0
    w = acos(dot(n,eVec)/norm(n)/ecc);
else
    w = 2*pi-acos(dot(n,eVec)/norm(n)/ecc);
end

%% True Anomaly
f = acos(dot(eVec,r0)/ecc/rMag);
if dot(r0,v0) <0
    f = 2*pi-f;
end