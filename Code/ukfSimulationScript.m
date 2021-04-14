clear all; clc; %close all
format longg
rng(123123);
% setup path
addpath('./Dynamics/')
addpath('./Measurements/')
addpath('./Utilities/')
addpath('./Filter/')
%% set parameters, settings, and initial conditions
Rbody = 250;%6378;

params.mu = 6.67430e-11*7.329e10;% m^3/s^2

% true initial conditions
a = 1.5*Rbody; ecc = 0; inc = 20*pi/180;w = 0; ra = 20*pi/180;f = 30*pi/180;
[r0,v0] = orbEl2rv(a,ecc,inc,w,ra,f,params.mu);
X0 = [r0;v0];

T = 2*pi*sqrt(a^3/params.mu);

params.L  = length(X0);
params.q  = 3;
params.r  = 3;
params.dt = 1;
params.options = odeset('RelTol',1e-8,'AbsTol',1e-8);
% landmark database (one landmark per column
params.numLand = 50;
params.landmark_db = Rbody*normc(random('Uniform',-1,1,params.r,params.numLand));
params.l = size(params.landmark_db,2);
params.P0 = diag([1000;1000;1000;1;1;1]);
params.procNoise = 1e-3*[10;10;10];% m^2/s^4
params.measNoise = 1/10*[1;1;1];% deg^2;

% ukf parameters
params.alpha     = 1e-4;
params.beta      = 2;
params.kappa     = 3 - params.L;
params.lambda    = params.alpha^2*(params.L+params.kappa)-params.L;
params.gamma     = sqrt(params.L+params.lambda);
params.W0mean    = params.lambda/(params.L + params.lambda);
params.W0cov     = params.W0mean + (1-params.alpha^2+params.beta);
params.Wi        = 1/(2*(params.L + params.lambda));
vProc = zeros(3,1);
params.dynamics = @(t,X) twoBodyEom(t,X,params.mu,vProc);

%% sensor simulation/trajectory simulation
% generate truth trajectory and measurements

% trajectory
t = 0:params.dt:T;
% measurements
% initialize arrays
losMeas = zeros(params.numLand,3,length(t));
Xtruth = zeros(length(t),params.L);
Xtruth(1,:) = X0(:)';
for ii = 2:length(t)
    
    %% trajectory
    vProc = mvnrnd(zeros(size(params.procNoise)),diag(params.procNoise));vProc =vProc(:);
    [~,Y] = ode45(@(t,X) twoBodyEom(t,X,params.mu,vProc),[t(ii-1),t(ii)],Xtruth(ii-1,:),params.options);
    Xtruth(ii,:) = Y(end,:);
    
    %% measurements
    trueMeas = meas_model_lod(Xtruth(ii,:)',1:params.numLand,params.landmark_db);
    
    for jj = 1:params.numLand
        % generate random rotation matrix to perturb measurement
        rpy = mvnrnd(zeros(3,1),diag(params.measNoise));% in degrees
        Rerror = Rz(rpy(3)*pi/180)*Ry(rpy(2)*pi/180)*Rx(rpy(1)*pi/180);
        
        losMeas(jj,:,ii) = Rerror*trueMeas(jj,:)';
    end
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
    fprintf('t=% f',t(ii))
    % UKF propagation step
    [Xfilter,Pfilter] = ukfPropStep(Xfilter,Pfilter,t(ii),t(ii-1),params);
    
    % loop through measurements
    for jj = 1:size(losMeas(:,:,ii),1)
        % UKF measurement step
        [yhatm,Pkyy,Pkxy] = ukfMeasStep(Xfilter,Pfilter,params,jj);
        
        % kalman gain
        K = Pkxy/Pkyy;
        
        % state update
        Xfilter = Xfilter + K*(losMeas(jj,:,ii)'-yhatm);
        
        % covariance update
        Pfilter = Pfilter - K*Pkyy*K';
    end
    
    % compute filter errors
    xEst(ii,:) = Xfilter;
    sig(ii,:) = sqrt(diag(Pfilter));
    fprintf('||posError = % f',norm(xEst(ii,1:3)-Xtruth(ii,1:3)));
    fprintf('||velError = % f',norm(xEst(ii,4:6)-Xtruth(ii,4:6)));
    fprintf('\n')
end
%% plotting and post sim processing
h = plotErrors(t,xEst,Xtruth,sig);