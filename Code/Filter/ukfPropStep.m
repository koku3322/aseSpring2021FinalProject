% Ken Kuppa, Dahlia Baker
% ASEN 6519
% Spring 2021
% last edited - KK, 4/11/2021

% performs one UKF propagation step
function [xhatm,Pkm] = ukfPropStep(X,P,tCurr,tPrev,dV,params)
dt = tCurr-tPrev;
%% IMU error compensation
dVCorr = dV - X(8:10);
%% propagation
% % make augmented covariance matrix
% Pkaug = zeros(params.L+length(params.imuBias.Sigma));
% Pkaug(1:params.L,1:params.L) = P;
% Pkaug(params.L+1:end,params.L+1:end) = params.imuBias.Sigma;
% get square root of the covariance matrix
Pksqrt = chol(P,'lower');

% generate sigma points using square root
kai = [X X+params.gamma*Pksqrt X-params.gamma*Pksqrt];

%%% propagate sigma points using non-linear dynamics
% initialize arrays
xhatm = zeros(params.L,1);
Pkm = zeros(params.L);
for ii = 1:2*params.L+1
    
    % for each sigma point, propagate state using compensated IMU Data
    r = norm(kai(1:3,ii));
    rhat = kai(1:3,ii)/r;
    g = -kai(7,ii)/r^2*rhat;
    vProc = g - dVCorr/dt;
    [~,Y] = ode45(@(t,X) twoBodyEom(t,X,vProc),[tPrev,tCurr],kai(:,ii),params.options);
    
    % assign final state to array
    kai(:,ii) = Y(end,:)';
    
    % compute mean state
    if ii == 1
        xhatm = xhatm + params.W0mean*kai(:,ii);
    else
        xhatm = xhatm + params.Wi*kai(:,ii);
    end
        
end
% reconstruct covariance using propagated sigma points
% covariance
for ii = 1:2*params.L+1
    
    % approximating a gaussian covariance using propagated sigma points and
    % computed mean
    if ii == 1
        Pkm = Pkm + params.W0cov*(kai(:,ii) - xhatm)*(kai(:,ii) - xhatm)';
    else
        Pkm = Pkm + params.Wi*(kai(:,ii) - xhatm)*(kai(:,ii) - xhatm)';
    end
    
end
% add process noise
Q = zeros(params.L);
Q(1:3,1:3) = (eye(3)+params.imuBias.Sigma)*dt;
Q(4:6,4:6) = (eye(3)+params.imuBias.Sigma);
Q(8:10,8:10) = (eye(3)+params.imuBias.Sigma);

Pkm = Pkm + Q;
Pkm = (Pkm + Pkm')/2;