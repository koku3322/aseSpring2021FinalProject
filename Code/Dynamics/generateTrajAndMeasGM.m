% generate truth trajectory and measurements
t = 0:params.dt:T;
% initialize arrays
losMeas = cell(length(t),2);
Xtruth = zeros(length(t),params.L);
Xtruth(1,:) = X0_true(:)';
for ii = 2:length(t)
    
    %% trajectory
    vProc = random(params.gmProc);vProc = vProc(:);
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
    
end
