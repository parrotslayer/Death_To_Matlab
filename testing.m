if UseReal == 1
        % Use real ECI with simulated errors
        Sat_ECI_est(:,t) = normrnd(Sat_ECI_true(:,t),sigma_orbit);
        
        % Convert ECI to ECEF
        Sat_ECEF_est = ECI_to_ECEF([Sat_ECI_est; t_since_equinox]);
    
        % Convert satellite's estimated ECEF position to est LLH
        Sat_LLH_est(t,:) = ECEF_to_LLH(Sat_ECEF_est(1),Sat_ECEF_est(2),Sat_ECEF_est(3));
    
    end
    %*************************************************************************
    
    % Set magnitude of ECI vector to Radius of the Earth
    Feature_ECI_est = Sat_ECI_est(:,t)*radius_earth/norm(Sat_ECI_est(:,t));
    
    % Convert ECI to ECEF
    Feature_ECEF_est = ECI_to_ECEF([Feature_ECI_est; t_since_equinox]);
       
    % Convert ECEF to LGCV wrt Satellite's position
    Feature_LGCV_est =  ECEF_to_LGCV(Sat_LLH_est(t,1),Sat_LLH_est(t,2),Sat_LLH_est(t,3),...
        Feature_ECEF_est(1), Feature_ECEF_est(2), Feature_ECEF_est(3));
    
    % Convert LGCV to Body using estimated angles
    Feature_Body_est = LGCV_to_Body(transpose(Attitude_est(t,:)),Feature_LGCV_est');
    
    % Add pointing errors to the Body frame coordinates
    pointing_errors = [0,0,0];
    pointing_errors = normrnd(pointing_errors,sigma_pointing_accuracy);
    Feature_Body_est = Rotate_3D(pointing_errors,Feature_Body_est);
    
    % Add ranging errors to the body readings
    Feature_Body_est = normrnd(Feature_Body_est,sigma_pointing_accuracy);
    
    % CALCULATE NADIR POINT
    % Convert Body to LGCV
    Nadir_LGCV_est = Body_to_LGCV(Attitude_est(t,:),Feature_Body_est);
    
    % Convert LGCV to ECEF
    Nadir_ECEF_est(t,:) = LGCV_to_ECEF(Sat_LLH_est(t,1),Sat_LLH_est(t,2),Sat_LLH_est(t,3),...
        Nadir_LGCV_est(1), Nadir_LGCV_est(2), Nadir_LGCV_est(3));
    
    % CALCULATE SWATH EDGE LOCATION
    % Rotate Body Frame by FOV/2 in the North Direction
    Swath_Body_est = C_rot*Feature_Body_est';
    
    % Scale Body Frame by D as distance is now different
    D = Viewing_Geometry(Sat_LLH_est(t,3),FOV/2);
    Swath_Body_est = Swath_Body_est/norm(Swath_Body_est)*D;
    
    % Body to LGCV
    Swath_LGCV_est = Body_to_LGCV(Attitude_est(t,:),Swath_Body_est');
    
    % LGCV to ECEF
    Swath_ECEF_est(t,:) =  LGCV_to_ECEF(Sat_LLH_est(t,1),Sat_LLH_est(t,2),Sat_LLH_est(t,3),...
        Swath_LGCV_est(1), Swath_LGCV_est(2), Swath_LGCV_est(3));
    
    % Calculate Nadir Errors
    Nadir_Error(t,:) = Nadir_ECEF_true(t,:) - Nadir_ECEF_est(t,:);
    Nadir_Error_Magnitude(t) = norm(Nadir_Error(t,:));
    
    % Calculate Swath Errors
    Swath_Error(t,:) = Swath_ECEF_true(t,:) - Swath_ECEF_est(t,:);
    Swath_Error_Magnitude(t) = norm(Swath_Error(t,:));
        
    % Get the Swath Width (difference between Nadir and Swath x2)
    Swath_Edge_est(t) = 2*norm(Swath_ECEF_est(t,:) - Nadir_ECEF_est(t,:));