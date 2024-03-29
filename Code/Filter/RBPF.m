%Rao-Blackwellized Particle Filter
%Ken Kuppa, Dahlia Baker
%ASEN 6519 Spring 2021
%Final Project
%last edited - DB 4/26/21

%note - much borrowed from UKF code

clear all; clc; close all
format longg
rng(123123);
% setup path
addpath('./Dynamics/')
addpath('./Measurements/')
addpath('./Utilities/')
addpath('./Filter/')
%%

%step 1 - setup
Rbody = 250;%6378;

params.mu = 6.67430e-11*7.329e10;% m^3/s^2

% true initial conditions
a = 1.2*Rbody; ecc = 0; inc = 20*pi/180;w = 0; ra = 20*pi/180;f = 30*pi/180;
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

%generate measurements and trajectory
vProc = zeros(3,1);
params.dynamics = @(t,X) twoBodyEom(t,X,params.mu,vProc);
generateTrajAndMeas();

%pf parameters
Qpf = diag([10,10,10,1e-2,1e-2,1e-2]);
R = params.procNoise;
Ns = 50;

%ukf parameters
params.alpha     = 1e-4;
params.beta      = 2;
params.kappa     = 3 - params.L;
params.lambda    = params.alpha^2*(params.L+params.kappa)-params.L;
params.gamma     = sqrt(params.L+params.lambda);
params.W0mean    = params.lambda/(params.L + params.lambda);
params.W0cov     = params.W0mean + (1-params.alpha^2+params.beta);
params.Wi        = 1/(2*(params.L + params.lambda));

%initialize particles and weights
xsamplehist = zeros(6,Ns,length(t));
Pksamplehist = zeros(6,6,Ns,length(t));
wsamplehist = (1/Ns)*ones(1,Ns);
weight_hist = zeros(1,Ns,length(t));
weight_hist(:,:,1) = wsamplehist;

xEst = nan(length(t),params.L,Ns);
sig = nan(length(t),params.L,Ns);
xEst(1,:,:) = X0.*ones(1,Ns);
%sig(1,:,:) = sqrt(diag(params.P0))';
updateApplied = zeros(size(t));

xsamplehist(:,:,1) = Xtruth(1,:)'.*ones(6,Ns,1);
for i = 1:Ns
    Pksamplehist(:,:,i,1) = params.P0;
    sig(1,:,i) = sqrt(diag(params.P0))';
end
%%
doplot = true;
%base filter - UKF
for i = 2:100%length(t)
   %for each particle
    for j = 1:Ns
       %prediction step
       %UKF
       %xfilter is the current particle, ran through a UKF
       Xfilter = xsamplehist(:,j,i-1)+mvnrnd(zeros(6,1),sqrt(Qpf))';
       Pfilter = Pksamplehist(:,:,j,i-1);
       [Xfilter,Pfilter] = ukfPropStep(Xfilter,Pfilter,t(i),t(i-1),params);
     
       %measurement update
       % pull out current measurement set from cell array
        losLand = losMeas{i,1};
        landIdx = losMeas{i,2};
        numObsLand = length(losMeas{i,2});
        % check if any measurement exist
        if ~isempty(losMeas{i,2})
            % loop through measurements
            for k = 1:numObsLand
                % UKF measurement step
                [yhatm,Pkyy,Pkxy] = ukfMeasStep(Xfilter,Pfilter,params,landIdx(k));
                % kalman gain
                K = Pkxy/Pkyy;
                % state update
                Xfilter = Xfilter + K*(losLand(:,k)-yhatm);
                % covariance update
                Pfilter = Pfilter - K*Pkyy*K';
                %recursive weight update
                wsamplehist(j) = sum(pdf('norm',losLand(:,k),yhatm,.1));
            end
            updateApplied(i) = true;
        else
            fprintf('No Meas')
        end
        %compute filter errors
        xEst(i,:,j) = Xfilter;
        sig(i,:,j) = sqrt(diag(Pfilter));
              
        xsamplehist(:,j,i) = Xfilter;
        Pksamplehist(:,:,j,i) = Pfilter;        
    end
    wsamplehist = wsamplehist./sum(wsamplehist);
    weight_hist(:,:,i) = wsamplehist;
    
    %Ness(i) = 1/(sum(weights(j,:).^2));
    %resampling step - put in later
    csw = cumsum(wsamplehist);
    
    urand = unifrnd(0,1,[1 Ns]);
    for m = 1:Ns
        %4 - move along the CSW
        l = 1;
        uj = urand(m);
        while uj > csw(l)
            l = l+1;
        end
        %5 - assign sample
        xsamplehist(:,m,i) = xsamplehist(:,l,i);
        %6 - assign weight
        %wsamplehist(m) = 1/Ns;
    end
    
    %plotting particles
    if doplot == true
        figure(1);
        hold on;
        plot(xsamplehist(1,:,i),xsamplehist(2,:,i),'ko');
        g = plot(Xtruth(i,1),Xtruth(i,2),'p','MarkerFaceColor','green','MarkerSize',15);
        xlim([100 600])
        ylim([200 700])
        legend([g],{'True position'});
        title('X-Y orbital plane: 50 particles','FontSize',16)
        xlabel('x (km)','FontSize',12)
        ylabel('y (km)','FontSize',12)
        hold off
        pause(0.00025)    
    else
        %do nothing
    end
end


h = plotErrors(t,xEst-Xtruth,sig,updateApplied);