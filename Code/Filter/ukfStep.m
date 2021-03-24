% Ken Kuppa, Dahlia Baker
% ASEN 6519
% Spring 2021
% last edited - KK, 3/23/2021

% performs one UKF step
function [xhatm,yhatm,Pkm,Pkyy,Pkxy] = ukfStep(X,P,tCurr,tPrev,params)

%% propagation
% get square root of the covariance matrix
Pksqrt = chol(P,'lower');

% generate sigma points using square root
kai = [X X+params.gamma*Pksqrt X-params.gamma*Pksqrt];

%%% propagate sigma points using non-linear dynamics
% initialize arrays
xhatm = zeros(params.n,1);
Pkm = zeros(params.n);
for ii = 1:2*params.L+1
    
    % for each sigma point, call ode45
    % NOTE: params.mu should be replaced by the estimated mu
    [~,Y] = ode45(dynamicsFunction,[tPrev,tCurr],kai(:,ii),params.mu);
    
    % assign final state to array
    kai(:,ii) = Y(:,end);
    
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
Pkm = Pkm + diag(params.procNoise);
%% update
% NOTE: could condition this section to only run one timesteps where we
% have measurements.

%%% generate sigma points for the predicted measurement
% initialize arrays
gammak = zeros(params.l,2*params.L+1);
yhatm = zeros(params.l,1);
Pkyy = zeros(params.l);
Pkxy = zeros(params.L,params.l);
for ii = 1:2*params.L+1
    
    % predicted measurement
    % NOTE: add landmark DB to params, ensure state indices are correct
    % add multiple landmark observation
    [los_vec] = meas_model_lod(kai(1:3,ii),landmarks_obs, params.landmark_db);
    
    
    gammak(:,ii) =  los_vec;
    % compute mean observation
    if ii == 1
        yhatm = yhatm + params.W0mean*gammak(:,ii);      
    else
        yhatm = yhatm + params.Wi*gammak(:,ii);
    end
    
end
% construct covariance using propagated state and observation sigma
for ii = 1:2*params.L+1
    if ii == 1
        Pkyy = Pkyy + params.W0cov*(gammak(:,ii)-yhatm)*(gammak(:,ii)-yhatm)';
        Pkxy = Pkxy + params.W0cov*(kaiProp(:,ii)-xhatm)*(gammak(:,ii)-yhatm)';
    else
        Pkyy = Pkyy + params.Wi*(gammak(:,ii)-yhatm)*(gammak(:,ii)-yhatm)';
        Pkxy = Pkxy + params.Wi*(kaiProp(:,ii)-xhatm)*(gammak(:,ii)-yhatm)';
    end
end
% add measurement noise
Pkyy = Pkyy + diag(params.measNoise);