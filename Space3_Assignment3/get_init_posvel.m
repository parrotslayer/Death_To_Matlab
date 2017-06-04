% Uses Herrick's Gibbs to estimate the velocity in ECI using atleast 3 ECI
% measurements.
% output = get_init_posvel(input,timeequinox,llhground)
% input = [time;range;azimuth;elevation] for multiple timesteps and are all
% read values (NaN filtered out previously)
% time_equinox (seconds)
% GS_LLH lat long height of ground station in degrees
% output = [r2;Vout] where r2 is the ECI coordinates of the 2nd data in
function output = get_init_posvel(input,time_equinox,GS_LLH)
%constants 
    G = 6.67408e-11 ; 
    M_earth = 5.972e24;
    u = G * M_earth;
    global time2
    
    %get vector of times for easy referencing
    time = input(1,:);  
    
    %RAE for first 3 observations
    obs1 = input(2:4,1);
    obs2 = input(2:4,2);
    obs3 = input(2:4,3);

    % PARAMETERS FOR HERRICK GIBBS
    % tij = tj - ti
    t23 = time(3) - time(2);
    t12 = time(2) - time(1);
    t13 = time(3) - time(1);
    
    %compute g1, g2, g3
    g1 = t23/(t12*t13);
    g3 = t12/(t23*t13);
    g2 = g1 - g3;
    
    %compute h1, h2, h3
    h1 = u * t23 / 12;
    h3 = u * t12 / 12;
    h2 = h1 - h3;
    
    % Compute the time since the equinox for each observation
    t_r1 = (time(1) - time_equinox);
    t_r2 = (time(2) - time_equinox);
    t_r3 = (time(3) - time_equinox);
    
    % Convert the polar LGCV observations to corresponding ECI coordinates
    % polar LGCV observations, at our ground station
    r1 = RAE_to_ECI(obs1,t_r1,GS_LLH); %remember to change the
    r2 = RAE_to_ECI(obs2,t_r2,GS_LLH); %time since earth rotates
    r3 = RAE_to_ECI(obs3,t_r3,GS_LLH); %slightly between obs.
    
    % Calculate intermediate vector variables for HG
    d1 = g1 + h1/(norm(r1))^3;
    d2 = g2 + h2/(norm(r2))^3;
    d3 = g3 + h3/(norm(r3))^3;
    
    % Get the velocity of the second observation
    Vout = - d1 * r1 + d2 * r2 + d3 * r3;
    
    output = [r2;Vout];
    time2 = time(2);
 end