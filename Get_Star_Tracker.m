% Get_Star_Tracker
% Takes in the satellite's position in LLH, the time since vernal equinox
% the star constellation in ECI and the Field of View in Radians
% Outputs the LGCV unit reference vector of all the visible satellites
% for the given tiemstep
%**********************************************************************
function output = Get_Star_Tracker(Sat_LLH,t_VE,Star_ECI,FOV )
% get number of stars
[num_stars, c] = size(Star_ECI);

%given star positions, need to convert to ECEF
Star_ECEF = zeros(num_stars,3);
for i = 1:num_stars
    Star_ECEF(i,:) = ECI_to_ECEF([Star_ECI(i,:),t_VE]);
end

% Number of visible satellites
vis = 0;
% Loop for all stars
for i = 1:num_stars
    %Convert ECEF position of the star to LGCV wrt the Satellite
    Star_LGCV = ECEF_to_LGCV(Sat_LLH(1),Sat_LLH(2),Sat_LLH(3)...
        , Star_ECEF(i,1), Star_ECEF(i,2), Star_ECEF(i,3));
    
    % Convert LGCV of Star to RAE (radians)
    Star_RAE = LGCV_to_RAE(Star_LGCV);
    
    % Check if elevation witin FOV
    if Star_RAE(3) > (pi/2 - FOV/2) && Star_RAE(3) < (pi/2 + FOV/2)
        %increment number of seen stars and store LGCV coordinates
        vis = vis + 1;
        %output normalised LGCV vector
        visible_stars(vis,:) = Star_LGCV/norm(Star_LGCV);
    end
end
output = visible_stars;
end
