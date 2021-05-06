% Ken Kuppa, Dahlia Baker
% ASEN 6519
% Spring 2021
% last edited - KK, 4/11/2021

% plots state errors with 2-sgima bounds. returns figure handle
function h = plotErrors(t,stateErr,sig,updateApplied)
h(1) = figure(1);
tiledlayout(3,2)
sgtitle({'State Estimation Errors:','Position & Velocity'})
for idx = 1:3
    % position
    ax(idx) = nexttile;
    plot(t,stateErr(:,idx),'b.-')
    hold on
    plot(t(~updateApplied),stateErr(~updateApplied,idx),'k.')
    plot(t,2*[-sig(:,idx),sig(:,idx)],'r')
    grid on, grid minor
    ylabel('Pos Error (m)')
    % velocity
    ax(idx+3) = nexttile;
    plot(t,stateErr(:,idx+3),'b.-')
    hold on
    plot(t(~updateApplied),stateErr(~updateApplied,idx+3),'k.')
    plot(t,2*[-sig(:,idx+3),sig(:,idx+3)],'r')
    grid on, grid minor
    ylabel('Vel Error (m/s)')
 
end
xlabel(ax(3),'Time (sec)')
xlabel(ax(6),'Time (sec)')
linkaxes(ax(1:3))
linkaxes(ax(4:6))

if size(stateErr,2)>6
    % plot mu errors
    h(2) = figure(2);
    plot(t,stateErr(:,7),'b.-')
    hold on
    plot(t(~updateApplied),stateErr(~updateApplied,7),'k.')
    plot(t,2*[-sig(:,7),sig(:,7)],'r')
    grid on, grid minor
    xlabel('Time (sec)')
    ylabel('Error (m^3/s^2)')
    title({'State Estimation Errors:','Gravitational Parameter (\mu)'})
    
    % imu errors
    h(3) = figure(3);
    plot(t,stateErr(:,8:10),'b.-')
    hold on
    plot(t(~updateApplied),stateErr(~updateApplied,8:10),'k.')
    plot(t,2*[-sig(:,8:10),sig(:,8:10)],'r')
    grid on, grid minor
    xlabel('Time (sec)')
    ylabel('Error (m/s)')
    title({'State Estimation Errors:','IMU Accel Bias'})
    
end