clear all; clc; close all
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
a = 2*Rbody; ecc = 0; inc = 20*pi/180;w = 0; ra = 20*pi/180;f = 30*pi/180;
[r0,v0] = orbEl2rv(a,ecc,inc,w,ra,f,params.mu);
X0 = [r0;v0];

T = 2*pi*sqrt(a^3/params.mu);

params.L  = length(X0);
params.q  = 3;
params.r  = 3;
params.dt = 5;
params.options = odeset('RelTol',1e-8,'AbsTol',1e-8);
% landmark database (one landmark per column)
params.numLand = 50;
azimuth = random('Uniform',0,pi,1,params.numLand);
elevation = random('Uniform',-pi/2,pi/2,1,params.numLand);
[x,y,z] = sph2cart(azimuth,elevation,Rbody);
params.landmark_db = [x;y;z];
params.P0 = diag([1000;1000;1000;1;1;1]);
params.procNoise = 1e-12*[10;10;10];% m^2/s^4
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
generateTrajAndMeas();
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
updateApplied = zeros(size(t));
for ii = 2:length(t)
    fprintf('t=% f',t(ii))
    % UKF propagation step
    [Xfilter,Pfilter] = ukfPropStep(Xfilter,Pfilter,t(ii),t(ii-1),params);
    
    % pull out current measurement set from cell array
    losLand = losMeas{ii,1};
    landIdx = losMeas{ii,2};
    numObsLand = length(losMeas{ii,2});
    
    % check if any measurement exist
    if ~isempty(losMeas{ii,2})
        % loop through measurements
        for jj = 1:numObsLand
            %
            % UKF measurement step
            [yhatm,Pkyy,Pkxy] = ukfMeasStep(Xfilter,Pfilter,params,landIdx(jj));
            
            % kalman gain
            K = Pkxy/Pkyy;
            
            % state update
            Xfilter = Xfilter + K*(losLand(:,jj)-yhatm);
            
            % covariance update
            Pfilter = Pfilter - K*Pkyy*K';
        end
        updateApplied(ii) = true;
    else
        fprintf('No Meas')
    end
    % compute filter errors
    xEst(ii,:) = Xfilter;
    sig(ii,:) = sqrt(diag(Pfilter));
    fprintf('||posError = % f',norm(xEst(ii,1:3)-Xtruth(ii,1:3)));
    fprintf('||velError = % f',norm(xEst(ii,4:6)-Xtruth(ii,4:6)));
    fprintf('\n')
end
%% plotting and post sim processing
h = plotErrors(t,xEst-Xtruth,sig,updateApplied);