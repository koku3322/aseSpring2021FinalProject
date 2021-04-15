%Dahlia Baker and Ken Kuppa
%ASEN 6519
%Spring 2021
%last edited - DB, 4/15/21

%ASE Final Project - Measurement Model
%Landmark Based OD and Gravity Estimation


%input - sc_state: the state of the spacecraft - [x,y,z,xdot,ydot,zdot]'.
%           We want these in the body-fixed frame, but that conversion
%           could be included in here later (it is not included now). We
%           don't really care about the velocity for now
%      - landmark_obs: this is a list of the landmarks observed, where the
%           vector entry corresponds to the database associated number
%      - landmark_db: the database of landmark positions in the
%           body-fixed frame l = [x,y,z,i] where i is the assigned number
%           of the landmark for db searching
%output - los_vec: stacked line of sight vectors, length i (number of 
%           observed landmarks) (difference in BF position
%           of s/c state and landmark in km)

function [los_vec] = meas_model_lod(sc_state, landmark_obs, landmark_db)

%pull observed landmarks from database
l_pos = landmark_db(1:3,landmark_obs);

% initialize array
los_vec = zeros(length(landmark_obs),3);
for i = 1:length(landmark_obs)
    los_vec(i,:) = sc_state(1:3)-l_pos(:,i);
end

% normalize 
los_vec = normr(los_vec);
end


