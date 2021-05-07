% generate truth trajectory and measurements
t = 0:params.dt:2*T;
% initialize arrays
losMeas = cell(length(t),2);
Xtruth = zeros(length(t),params.L);
Xtruth(1,:) = X0_true(:)';
deltaVTrue = zeros(3,length(t));
deltaV = zeros(3,length(t));
X = zeros(size(Xtruth));
X(1,:) = X0_true(:)'; 
for ii = 2:length(t)
    
    %% trajectory
    vProc = zeros(3,1);
    [~,Y] = ode45(@(t,X) twoBodyEom(t,X,vProc),[t(ii-1),t(ii)],Xtruth(ii-1,:),params.options);
    Xtruth(ii,:) = Y(end,:);
    
    %% measurements
    trueMeas = meas_model_lod(Xtruth(ii,:)',1:params.numLand,params.landmark_db);
    % check which landmarks are visible
    visIdx = find((dot(normc(params.landmark_db),trueMeas'))>0);
        
    if ~isempty(visIdx)
        currMeas = zeros(3,length(visIdx));
        for jj = 1:length(visIdx)
            % generate random rotation matrix to perturb measurement
            rpy = random(params.gmMeas);
%             rpy = mvnrnd(zeros(3,1),diag(params.measNoise));% in degrees
            Rerror = Rz(rpy(3)*pi/180)*Ry(rpy(2)*pi/180)*Rx(rpy(1)*pi/180);
            
            currMeas(:,jj) = Rerror*trueMeas(visIdx(jj),:)';
        end
    else
        currMeas = [];
    end
    % assign set of measurements to cell array
    losMeas{ii,1} = currMeas;
    
    % assign index of observed landmarks
    losMeas{ii,2} = visIdx;
    
    
    %% IMU data
    
    % truth accel data
    deltaVTrue(:,ii) = Xtruth(ii,4:6)-Xtruth(ii-1,4:6);% actual difference in velocity
    
%     % add error
    biasErr = random(params.imuBias);biasErr = biasErr(:);

    deltaV(:,ii) = deltaVTrue(:,ii) + biasErr;
        
    vel = X(ii-1,4:6) + deltaV(:,ii)';
    pos = X(ii-1,1:3) + vel*params.dt + 0.5*deltaV(:,ii)'*params.dt;
    
    
    X(ii,1:7) = [pos,vel,params.mu];
    X(ii,8:10) = sum(params.imuBias.mu.*wProc)';
end
