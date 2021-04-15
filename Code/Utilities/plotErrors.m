% Ken Kuppa, Dahlia Baker
% ASEN 6519
% Spring 2021
% last edited - KK, 4/11/2021

% plots state errors with 2-sgima bounds. returns figure handle
function h = plotErrors(t,stateErr,sig,updateApplied)
h = figure;
tiledlayout(3,2)
sgtitle('State Estimation Errors')
for idx = 1:3
    % position
    ax(idx) = nexttile;
    plot(t,stateErr(:,idx),'k.-')
    hold on
    plot(t(~updateApplied),stateErr(~updateApplied,idx),'r.')
    plot(t,2*[-sig(:,idx),sig(:,idx)],'r')
    grid on, grid minor
    ylabel('Pos Error (m)')
    % velocity
    ax(idx+3) = nexttile;
    plot(t,stateErr(:,idx+3),'k.-')
    hold on
    plot(t(~updateApplied),stateErr(~updateApplied,idx+3),'r.')
    plot(t,2*[-sig(:,idx+3),sig(:,idx+3)],'r')
    grid on, grid minor
    ylabel('Vel Error (m/s)')
 
end
xlabel(ax(3),'Time (sec)')
xlabel(ax(6),'Time (sec)')
linkaxes(ax(1:3))
linkaxes(ax(4:6))