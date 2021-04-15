% Ken Kuppa, Dahlia Baker
% ASEN 6519
% Spring 2021
% last edited - KK, 4/14/2021

% performs one UKF measurement step
function [yhatm,Pkyy,Pkxy] = ukfMeasStep(xhatm,Pkm,params,landIdx)


%%% generate sigma points for the predicted measurement
Pksqrt = chol(Pkm,'lower');
% generate sigma points using square root
kai = [xhatm xhatm+params.gamma*Pksqrt xhatm-params.gamma*Pksqrt];

% initialize arrays
gammak = zeros(params.r,2*params.L+1);
yhatm = zeros(params.r,1);
Pkyy = zeros(params.r);
Pkxy = zeros(params.L,params.r);
for ii = 1:2*params.L+1
    
    % predicted measurement
    % NOTE: add landmark DB to params, ensure state indices are correct
    % add multiple landmark observation
    [gammak(:,ii)] = meas_model_lod(kai(1:3,ii),landIdx, params.landmark_db);
    
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
        Pkxy = Pkxy + params.W0cov*(kai(:,ii)-xhatm)*(gammak(:,ii)-yhatm)';
    else
        Pkyy = Pkyy + params.Wi*(gammak(:,ii)-yhatm)*(gammak(:,ii)-yhatm)';
        Pkxy = Pkxy + params.Wi*(kai(:,ii)-xhatm)*(gammak(:,ii)-yhatm)';
    end
end
% add measurement noise
R = Rz(params.measNoise(3)*pi/180)*Ry(params.measNoise(2)*pi/180)*Rx(params.measNoise(1)*pi/180);
Pkyy = Pkyy + R;