clear; close all; clc
format longg
rng(123123);
% setup path
addpath('./Dynamics/')
addpath('./Measurements/')
addpath('./Utilities/')
addpath('./Filter/')

%% setup initial condition
a = 6800; ecc = 0; inc = 90*pi/180;w = 0; ra = 20*pi/180;f = 30*pi/180; mu = 398600;
[r0,v0] = orbEl2rv(a,ecc,inc,w,ra,f,mu);
state0 = [r0;v0];

T = 10*pi*sqrt(a^3/mu);
% tSpan = linspace(0,10*T,10000);
tSpan = 0:.1:T;
%% Plot earth
Rearth = 6378;
[xE,yE,zE] = sphere(50);
surf(Rearth*xE,Rearth*yE,Rearth*zE,'EdgeColor','none','FaceColor','b');
hold on; grid on, grid minor
axis equal

%% propagate using point mass dynamics
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[~,YPoint] = ode45(@(t,X) twoBodyEom(t,X,mu),tSpan,state0,options);
plot3(YPoint(:,1),YPoint(:,2),YPoint(:,3),'r')


%% propagate using oblateness dynamics
J2 = 1.082636e-3;
[~,YJ2] = ode45(@(t,X) twoBodyJ2Eom(t,X,mu,J2,Rearth,0),tSpan,state0,options);
plot3(YJ2(:,1),YJ2(:,2),YJ2(:,3),'g')