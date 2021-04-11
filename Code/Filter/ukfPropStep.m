% Ken Kuppa, Dahlia Baker
% ASEN 6519
% Spring 2021
% last edited - KK, 4/11/2021

% performs one UKF propagation step
function [xhatm,Pkm] = ukfPropStep(X,P,tCurr,tPrev,params)

%% propagation
% get square root of the covariance matrix
Pksqrt = chol(P,'lower');

% generate sigma points using square root
kai = [X X+params.gamma*Pksqrt X-params.gamma*Pksqrt];

%%% propagate sigma points using non-linear dynamics
% initialize arrays
xhatm = zeros(params.L,1);
Pkm = zeros(params.L);
for ii = 1:2*params.L+1
    
    % for each sigma point, call ode45
    % NOTE: params.mu should be replaced by the estimated mu
    [~,Y] = ode45(params.dynamics,[tPrev,tCurr],kai(:,ii),params.mu);
    
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
Pkm = Pkm + diag(params.procNoise);