 function params = nlls_orbit_determ(vernal_s,obs,GS_ECEF,init_posvel_guess,gs_llh)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NLLS_ORBIT_DETERM - Uses Non-linear Least Squares technique to estimate
% the orbital parameters given a series of ground station tracking
% obserations.
%
% Inputs: obs - N x 4 vector of [time, range, azimuth, elevation] measurements
%         from a single ground station
%         GS_ECEF - ECEF position of the ground tracking station
%         init_posvel_guess - Initial parameter estimate 
%         vernal_s - The date at which TLE code taken in seconds
%         gs_llh - LLH coordinates of the ground stations
% Outputs: params - Orbital parameters: [a,e,i,RAAN,AoP,Mo(at epoch)]
%
% NOTE: Elements of this function are missing! you will need to fill in
% these gaps to get the code working!

% This version of the code is designed for a single ground station; you
% will need to make several adjustments for multiple ground stations 

% See Week 6 Slides 19 to 31 for Details on non-linear least squares for
% orbit determination from ground station tracking measurements
% [10951792.6290923;-20903284.0364744;-8160479.16933671;2111.69338749256;-211.780587177724;3434.47872583941]
% FILL IN HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialise H matrix and error vector
% Calculate the total number of ground station observations
num_obs = size(obs,1); 
% Initialise the total measurement Jacobian matrix
H = zeros(3*num_obs,6); 
% Initialise the residual vector
delta_y = zeros(3*num_obs,1); 
% Use initial guess at orbital parameters to get initial position and
% velocity 
X = init_posvel_guess;
delta_x = 10e9;
% Main iteration loop to converge on initial satellite position and
% velocity
max_iter = 20;
tol = 1e-6;
iter = 0;
while (norm(delta_x) > tol)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Re-linearise the system, build the observation model and Jacobian
    
    % Loop through each observation, build observation matrix
    for i = 1:num_obs
        
        % Use current value of predicted initial position to calculate
        % position at time of observation
        ini_time = obs(2,1);
        current_obs_time = obs(i,1);
        % current_gs = obs(i,5);
        delta_t = current_obs_time - ini_time;
        [pos, vel] = universal_conic_section_orbit(delta_t, X(1:3), X(4:6));
        % Get the expected measurement (range, azimuth, elevation)
        % FILL IN HERE: expected measurement as a function of ground
        % station location, and satellite location
        time_VE = current_obs_time - vernal_s;  %time since VE passage
        %ECI -> ECEF
        ECEF = ECI_to_ECEF([pos;time_VE]);
        % ECEF + LLH -> LGCV
        LGCV = ECEF_to_LGCV(gs_llh(1), gs_llh(2), gs_llh(3), ECEF(1), ECEF(2), ECEF(3));
        %LGCV -> RAE
        expected_obs = LGCV_to_RAE(LGCV);

        % Build H matrix
        % Time since the first observation
        dxdxo = universal_conic_section_jacobian(pos, vel, X(1:3), X(4:6),current_obs_time -  ini_time);
        dhdx = rangeazielev_obs_jacob(GS_ECEF, pos, time_VE);
        H((3*i - 2):3*i,:) = dhdx*dxdxo; % Builds the rows of the Jacobian matrix corresponding to this measurement
        
        % Build the residual vector for this measurement
        actual_obs = obs(i,2:4);
        delta_y((3*i - 2):3*i,1) = (actual_obs - expected_obs')';
    end
    
    % Use non-linear least squares to estimate error in x
 %    delta_x = inv(H'*H)*H'*delta_y;
    
    % Alternative formulation (much better, type "help mldivide" for
    % details)
    delta_x = (H'*H)\H'*delta_y;
    
    X = X + delta_x;
    
    iter = iter + 1;
    
    if (iter >= max_iter)
        disp('Failed to Converge !!')
        break;
    end    
end

% FILL IN HERE: Transform iterated estimate of initial position and
% velocity into orbital parameters (See Notes Week 3, Week 6 Slide 31)
% use Herrick's Gibbs to convert Vel and Pos into orbital parameters
params = ECI_posvel_to_Orbit(X);
