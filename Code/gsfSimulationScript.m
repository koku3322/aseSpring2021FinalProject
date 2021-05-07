clear all; clc; close all
format longg
rng(123123);
% setup path
addpath('./Dynamics/')
addpath('./Measurements/')
addpath('./Utilities/')
addpath('./Filter/')
gcp;
%% set parameters, settings, and initial conditions
Rbody = 250;%6378e3;

params.mu = 6.67430e-11*7.329e10;% m^3/s^2

% imu errors
wProc = 1/3*ones(3,1);
mean1  = 1e-8*[1,1,1];
mean2  = 2e-8*[1,1,1];
mean3  = -1e-8*[1,1,1];
sigma = diag(1e-10*[1;1;1]);
params.imuBias = gmdistribution([mean1;mean2;mean3],sigma,wProc);

% true initial conditions (5 mixands)
a = 1.2*Rbody; ecc = 0; inc = 20*pi/180;w = 0; ra = 20*pi/180;f = 20*pi/180;
[r0,v0] = orbEl2rv(a,ecc,inc,w,ra,f,params.mu);
X0_true1 = [r0;v0;4.5;params.imuBias.mu(1,:)'];

a = 1.2*Rbody; ecc = 0; inc = 20*pi/180;w = 0; ra = 20*pi/180;f = 25*pi/180;
[r0,v0] = orbEl2rv(a,ecc,inc,w,ra,f,params.mu);
X0_true2 = [r0;v0;5;params.imuBias.mu(2,:)'];

a = 1.2*Rbody; ecc = 0; inc = 20*pi/180;w = 0; ra = 20*pi/180;f = 30*pi/180;
[r0,v0] = orbEl2rv(a,ecc,inc,w,ra,f,params.mu);
X0_true3 = [r0;v0;5.5;params.imuBias.mu(3,:)'];

a = 1.2*Rbody; ecc = 0; inc = 20*pi/180;w = 0; ra = 20*pi/180;f = 35*pi/180;
[r0,v0] = orbEl2rv(a,ecc,inc,w,ra,f,params.mu);
X0_true4 = [r0;v0;4;params.imuBias.mu(3,:)'];

a = 1.2*Rbody; ecc = 0; inc = 20*pi/180;w = 0; ra = 20*pi/180;f = 40*pi/180;
[r0,v0] = orbEl2rv(a,ecc,inc,w,ra,f,params.mu);
X0_true5 = [r0;v0;4;params.imuBias.mu(3,:)'];

T = 2*pi*sqrt(a^3/params.mu);

params.L  = length(X0_true1);
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
params.P0 = zeros(params.L);
params.P0(1:7,1:7) = diag([1000;1000;1000;1;1;1;1]);
params.P0(8:10,8:10) = params.imuBias.Sigma;
params.measNoise = [.1;.1;.1];% deg^2;

% make gm prior
wPrior = 1/3*ones(5,1);%[0.1;0.8;0.1];
gmPrior = gmdistribution([X0_true1,X0_true2,X0_true3,X0_true4,X0_true5]',params.P0,wPrior);

% make gm meas noise
wMeas = 1/3*ones(3,1);
measNoise(:,:,1) = diag([0.1;0.1;0.1]);
measNoise(:,:,2) = diag([1;1;1]);
measNoise(:,:,3) = diag([0.5;0.5;0.5]);
params.gmMeas = gmdistribution(zeros(3,size(measNoise,3)),measNoise,wMeas);

% mixture compression
params.Mdes = 5;
params.compressionMethod = 'brutalTrunc';%'momentMatch';%
params.nCompSkip = 1;

% ukf parameters
params.alpha     = 1e-4;
params.beta      = 2;
params.kappa     = 3 - params.L;
params.lambda    = params.alpha^2*(params.L+params.kappa)-params.L;
params.gamma     = sqrt(params.L+params.lambda);
params.W0mean    = params.lambda/(params.L + params.lambda);
params.W0cov     = params.W0mean + (1-params.alpha^2+params.beta);
params.Wi        = 1/(2*(params.L + params.lambda));
%% sensor simulation/trajectory simulation
X0_true = X0_true1;%random(gmPrior);
generateTrajAndMeasGM();
%% filter
% initialize filter
% Xprior = gmPrior.mu';
% Pprior = repmat(params.P0,1,1,3); 
Xprior = cell(gmPrior.NumComponents,1);
Pprior = cell(gmPrior.NumComponents,1);
for ii = 1:length(Xprior)
    Xprior{ii,:} = gmPrior.mu(ii,:)';
    Pprior{ii,:} = params.P0;
end
% initialize arrays
xEst = nan(length(t),params.L);
sig = nan(length(t),params.L);
xEst(1,:) = sum(gmPrior.mu'.*wPrior',2);
sig(1,:) = sqrt(diag(params.P0))';
updateApplied = zeros(size(t));
% params.P0 = params.P0*10;
% params.procNoise = 1e-8*[10;10;10];% m^2/s^4
% params.measNoise = 2*[1;1;1];% deg^2;
wts = wPrior;
for ii = 2:length(t)
    fprintf('t=% f',t(ii))
    %% GSUKF propagation step
    [Xprior,Pprior] = gsukfPropStep(Xprior,Pprior,t(ii),t(ii-1),deltaV(:,ii),params);
    
    % pull out current measurement set from cell array
    losLand = losMeas{ii,1};
    landIdx = losMeas{ii,2};
    numObsLand = length(losMeas{ii,2});
    
    % update weights for the propagated mixands
    wts = wts*wProc';wts = wts(:);
    
    %% check if any measurement exist
    if numObsLand ~= 0
        muProp = repmat(Xprior,1,length(params.gmMeas.mu));
        Pprop  = repmat(Pprior,1,length(params.gmMeas.mu));
        % updated weights
        wtsUpd = repmat(wts,1,length(params.gmMeas.mu));
        
        % GSUKF update step
        for iPrior = 1:size(muProp,1)
            muPrior_i = Xprior{iPrior};
            Pprior_i = Pprior{iPrior};
            
            parfor iMeas = 1:size(muProp,2)
                
                % loop through measurements
                for kk = 1:numObsLand

                    % UKF measurement step
                    [yhatm,Pkyy,Pkxy] = gsukfMeasStep(muProp{iPrior,iMeas},Pprop{iPrior,iMeas},params,landIdx(kk));
                    
                    measNoise = params.gmMeas.Sigma(:,:,iMeas);
                    R = Rz(measNoise(3)*pi/180)*Ry(measNoise(2)*pi/180)*Rx(measNoise(1)*pi/180);
                    Pkyy = Pkyy + R;
                    Pkyy = (Pkyy+Pkyy')/2;
                    % skip measurement if non PSD
                    if min(eig(Pkyy)) < 0
%                         fprintf('skipping measurement % i',kk);
                        break
                    end
                    % kalman gain
                    K = Pkxy/Pkyy;

                    % state update
                    muProp{iPrior,iMeas} = muProp{iPrior,iMeas} + K*(losLand(:,kk)-yhatm);

                    % covariance update
                    P = Pprop{iPrior,iMeas} - K*Pkyy*K';
                    Pprop{iPrior,iMeas} = (P + P')/2;

                    % update weights
                    wtsUpd(iPrior,iMeas) = ...
                        wtsUpd(iPrior,iMeas)...
                        ...*params.gmMeas.ComponentProportion(iMeas)...
                        *mvnpdf(losLand(:,kk),yhatm,Pkyy);

                end
            end
            
        end
        updateApplied(ii) = true;
        
        % make output into colum cell vector
        Xprior = muProp(:);
        Pprior  = Pprop(:);
        
        % normalize the weights
        wtsUpd = wtsUpd(:);
        wts = wtsUpd/sum(wtsUpd);
        
        
    end
    
    
    %% compute MMSE error and covariance
    xhat = zeros(params.L,1);
    Phat = zeros(params.L);
    for jj = 1:length(Xprior)
        % mean
        xhat = xhat + wts(jj).*Xprior{jj};
        
        % covariance
        Pbar = Pprior{jj} + Xprior{jj}*Xprior{jj}';
        Phat = Phat + wts(jj)*Pbar;
        
    end
    Phat = Phat - xhat*xhat';
    Phat = (Phat+Phat')/2;
    
    %% mixture compression
    
    switch params.compressionMethod
        case 'momentMatch'% moment matching
            Xprior = {xhat};
            Pprior = {Phat};
            wts = 1;
        case 'brutalTrunc'% brutal truncation
            if mod(ii,params.nCompSkip)==0
                [wts,Xprior,Pprior] = brutalTrunc(wts,Xprior,Pprior,params.Mdes);
            end
    end
    %% compute filter errors
    xEst(ii,:) = xhat(:)';
    sig(ii,:) = sqrt(diag(Phat));
    fprintf('||numObsLand = % i',numObsLand);
    fprintf('||posError = % f',norm(xEst(ii,1:3)-Xtruth(ii,1:3)));
    fprintf('||velError = % f',norm(xEst(ii,4:6)-Xtruth(ii,4:6)));
    fprintf('\n')
end
%% plotting and post sim processing
switch params.compressionMethod
    case 'momentMatch'% moment matching
        tag = [params.compressionMethod];
    case 'brutalTrunc'% brutal truncation
        tag = [params.compressionMethod,'_mdes',num2str(params.Mdes)];
end
h = plotErrors(t,xEst-Xtruth,sig,updateApplied);
saveas(h(1),['.\..\Results\posVelErrs_gsf_',tag],'png')
saveas(h(1),['.\..\Results\posVelErrs_gsf_',tag],'fig')
if length(h)>1
    saveas(h(2),['.\..\Results\muErrs_gsf_',tag],'png')
    saveas(h(2),['.\..\Results\muErrs_gsf_',tag],'fig')
    saveas(h(3),['.\..\Results\imuErrs_gsf_',tag],'png')
    saveas(h(3),['.\..\Results\imuErrs_gsf_',tag],'fig')
end

save(['.\..\Results\filterdata_',tag],'t','xEst','sig')