%% Get True Position of the Nadir Point and North Swath Edge
% Convert satellite's estimated ECEF position to est LLH
Sat_LLH_true(t,:) = ECEF_to_LLH(Sat_ECEF_true(1,t),Sat_ECEF_true(2,t),Sat_ECEF_true(3,t));

% Set magnitude of ECI vector to Radius of the Earth
Feature_ECI_true = Sat_ECI_true(:,t)*radius_earth/norm(Sat_ECI_true(:,t));

% Convert ECI to ECEF
current_time = time_epoch + t;    %time since last epoch for time = n
t_since_equinox = current_time - julian_date17;
Nadir_ECEF_true(t,:) = ECI_to_ECEF([Feature_ECI_true; t_since_equinox]);

% Convert ECEF to LGCV
Feature_LGCV_true =  ECEF_to_LGCV(Sat_LLH_true(t,1),Sat_LLH_true(t,2),Sat_LLH_true(t,3),...
    Nadir_ECEF_true(t,1), Nadir_ECEF_true(t,2), Nadir_ECEF_true(t,3));  

% Convert LGCV to Body
Feature_Body_true = LGCV_to_Body(Attitude_Real(:,t),Feature_LGCV_true');

% Rotate Body Frame by FOV/2 in the North Direction
th = FOV/2;
C_rot = [cos(th) 0 -sin(th);
    0 1 0;
    sin(th) 0 cos(th)];
Feature_Body_true = C_rot*Feature_Body_true';

% Body to LGCV
Feature_LGCV_true = Body_to_LGCV(Attitude_Real(:,t),Feature_Body_true');

% LGCV to ECEF
temp =  LGCV_to_ECEF(Sat_LLH_true(t,1),Sat_LLH_true(t,2),Sat_LLH_true(t,3),...
    Feature_LGCV_true(1), Feature_LGCV_true(2), Feature_LGCV_true(3)); 

