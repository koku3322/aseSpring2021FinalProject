clear all; clc; close all
format longg
rng(123123);
% setup path
addpath('./Dynamics/')
addpath('./Measurements/')
addpath('./Utilities/')
addpath('./Filter/')
%% set parameters, settings, and initial conditions


params.L = 6;
params.q = 3;
params.r = 3;
params.mu        = 3.986e5;% km^3/s^2
params.dt = 1;
params.options = odeset('RelTol',1e-8,'AbsTol',1e-8);
% landmark database (one landmark per column
params.landmark_db = [6378,0,0,0,6378*cosd(20);
                      0,6378,0,6378*cosd(45),0;
                      0,0,6378,6378*sind(45),6378*sind(20)];% km
params.numLand = size(params.landmark_db,2);
params.l = size(params.landmark_db,2);
params.P0 = diag([1;1;1;.01;.01;.01]);
params.procNoise = [.1;.1;.1;.001;.001;.001];
params.measNoise = [10;10;10];

% ukf parameters
params.alpha     = 1e-4;
params.beta      = 2;
params.kappa     = 3 - params.L;
params.lambda    = params.alpha^2*(params.L+params.kappa)-params.L;
params.gamma     = sqrt(params.L+params.lambda);
params.W0mean    = params.lambda/(params.L + params.lambda);
params.W0cov     = params.W0mean + (1-params.alpha^2+params.beta);
params.Wi        = 1/(2*(params.L + params.lambda));
params.dynamics = @(t,X) twoBodyEom(t,X,params.mu);

% true initial conditions
a = 6800; ecc = 0; inc = 20*pi/180;w = 0; ra = 20*pi/180;f = 30*pi/180;
[r0,v0] = orbEl2rv(a,ecc,inc,w,ra,f,params.mu);
X0 = [r0;v0];

T = 2*pi*sqrt(a^3/params.mu);
%% sensor simulation/trajectory simulation
% generate truth trajectory and measurements

% trajectory
t = 0:params.dt:T;
[~,Xtruth] = ode45(params.dynamics,t,X0,params.options);
% add process noise
Xtruth = Xtruth + mvnrnd(zeros(6,1),diag(params.procNoise),length(t));
% measurements
trueLos = zeros(params.numLand,3,length(t));
for ii = 1:length(t)
    
    trueLos(:,:,ii) = meas_model_lod(Xtruth(ii,:)',1:params.numLand,params.landmark_db)+...
                        mvnrnd(zeros(3,1),diag(params.measNoise));
    trueLos(:,:,ii) = normr(trueLos(:,:,ii));
end

%% filter
% initialize filter
X0 = mvnrnd(Xtruth(1,:),params.P0);
Xfilter = X0(:);
Pfilter = params.P0;
% initialize arrays
xEst = nan(length(t),params.L);
sig = nan(length(t),params.L);
xEst(1,:) = X0;
sig(1,:) = sqrt(diag(params.P0))';
for ii = 2:length(t)
    
    % UKF propagation step
    [Xfilter,Pfilter] = ukfPropStep(Xfilter,Pfilter,t(ii),t(ii-1),params);
    
    % loop through measurements
    for jj = 1:size(trueLos(:,:,ii),1)
        % UKF measurement step
        [yhatm,Pkyy,Pkxy] = ukfMeasStep(Xfilter,Pfilter,params,jj);
        
        % kalman gain
        K = Pkxy/Pkyy;
        
        % state update
        Xfilter = Xfilter + K*(trueLos(jj,:,ii)'-yhatm);
        
        % covariance update
        Pfilter = Pfilter - K*Pkyy*K';
    end
    
    % compute filter errors
    xEst(ii,:) = Xfilter;
    sig(ii,:) = sqrt(diag(Pfilter));
end
%% plotting and post sim processing
close all
h = plotErrors(t,xEst,Xtruth,sig);