% Ken Kuppa, Dahlia Baker
% ASEN 6519
% Spring 2021
% last edited - KK, 4/11/2021

% performs mixture compression using brutal truncation
function [wts,X,P] = brutalTrunc(wts,X,P,Mdes)
if length(X)<=Mdes
    fprintf('skipping compression\n');
    return
end
[wts,idx]=maxk(wts,Mdes);

% renormalize weights
wts = wts/sum(wts);

% update prior and covariance
X = X(idx);
P = P(idx);

