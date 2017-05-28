%% Get Estimated Nadir Point
% Set magnitude of ECI vector to Radius of the Earth
Feature_ECI_est = Sat_ECI_est(:,t)*radius_earth/norm(Sat_ECI_est(:,t));

% Convert ECI to ECEF
Feature_ECEF_est = ECI_to_ECEF([Feature_ECI_est; t_since_equinox]);

% Convert ECEF to LGCV wrt Satellite's position
Feature_LGCV_est =  ECEF_to_LGCV(Sat_LLH_est(t,1),Sat_LLH_est(t,2),Sat_LLH_est(t,3),...
    Feature_ECEF_est(1), Feature_ECEF_est(2), Feature_ECEF_est(3));  

% Convert LGCV to Body using estimated angles
Feature_Body_est = LGCV_to_Body(transpose(Attitude_Est(t,:)),Feature_LGCV_est');

% Add Errors to Body readings
Nadir_Body_est = normrnd(Feature_Body_est,sigma_tracking_errors);

% Convert Body to LGCV
Nadir_LGCV_est = Body_to_LGCV(Attitude_Est(t,:),Nadir_Body_est);

% Convert LGCV to ECEF
Nadir_ECEF_est(t,:) = LGCV_to_ECEF(Sat_LLH_est(t,1),Sat_LLH_est(t,2),Sat_LLH_est(t,3),...
    Nadir_LGCV_est, Nadir_LGCV_est, Nadir_LGCV_est);

% Calculate Nadir Errors
Nadir_Error(t,:) = Nadir_ECEF_true(t,:) - Nadir_ECEF_est(t,:);
Nadir_Error_Magnitude(t) = norm(Nadir_Error(t,:));