% AERO4701 Space 3 Assignment 2 Question 1
clc
%close all
clear all

%% Part A
%Constants
rad = 6378137; %earth radius km
deg2rad = pi/180;
rad2deg = 180/pi;
G = 6.67408e-11;
M_earth = 5.972e24;
mu = G*M_earth;
%R_ECI = zeroes(size(R_ECI_2D, 1), size(R_ECI_2D, 2), 0);

% Import Data
GPS_data = dlmread('GPSsat_ephem.txt');

%loop for each satellite
for s = 1:31
    
    %% Extract orbital parameters
    a = GPS_data(s,2);          %m
    e = GPS_data(s,3);
    inc = GPS_data(s,4)*deg2rad;        % (rads)
    Omega = GPS_data(s,5)*deg2rad;      %Right ascension of ascending node (rads)
    omega = GPS_data(s,6)*deg2rad;      %Argument of Perigee (rads)
    M = GPS_data(s,7)*deg2rad;          %Mean anomaly at epoch (rads)
    t_epoch = GPS_data(s,8);    %seconds
    
    
    %% calculate Orbital elements
    n = sqrt(mu/a^3);
    period = 86400/n;    % seconds
    p = a*(1 - e^2);
    tau = 1 - M/n;
    
    t_VE2016 = [2016 09 23 12 21 0];    %date of last VE
    t_VEepoch = 7347737.366;            %time since last VE from assignment
    t_VE = t_epoch - t_VEepoch;         %time elapsed between TLE and last VE
    
    %% Get ECI for time steps
    
    %looping parameters go here
    i = 1;
    t = 0;
    dt = 100;
    
    for t=0:dt:43200 
        %% Solve initial orbital parameters for t
        %solve for M
        M = n*(t - tau);
        
        % solver keplers equation to numerically find E
        tol = 0.001;
        E = solveE(M,e,tol);
        
        %solve for theta/true anomaly
        theta = 2*atan(((1+e)/(1-e))^0.5*tan(E/2));
        
        %solve for radius
        r = p/(1+e*cos(theta));
        
        %% Get ECI,ECEF,LLH positions for the time step
        % Transform orbit coordinate to ECI, Variable "R_Orbit"
        xi = r.*cos(theta);
        yi = r.*sin(theta);
        zi = 0;
        
        % Coordinate Transformation Matrix
        Orbit2ECI = [(cos(Omega)*cos(omega)-sin(Omega)*sin(omega)*cos(inc)) (-cos(Omega)*sin(omega)-sin(Omega)*cos(omega)*cos(inc)) (sin(Omega)*sin(inc));
            (sin(Omega)*cos(omega)+cos(Omega)*sin(omega)*cos(inc)) (-sin(Omega)*sin(omega)+cos(Omega)*cos(omega)*cos(inc)) (-cos(Omega)*sin(inc));
            (sin(omega)*sin(inc)) cos(omega)*sin(inc) cos(inc)];
        
        % Create get the ECI orbit for the current time step
        % R_ECI = ECI_C_Orbit * R_Orbit
        temp = Orbit2ECI*[xi;yi;zi];
        
        %copy results to R_ECI matrix
        R_ECI(1,i,s) = temp(1);   %X
        R_ECI(2,i,s) = temp(2);   %Y
        R_ECI(3,i,s) = temp(3);   %Z
        
        %increment looping variables
        i = i + 1;
        
    end
end


%% Plot 3D ECI Orbit

figure
load('topo.mat','topo','topomap1');

% Create a sphere, make it earth sized (in meters)
[x,y,z] = sphere(50);
x = x.*6378000;
y = y.*6378000;
z = z.*6378000;

% Set visual properties for plot
props.AmbientStrength = 0.1;
props.DiffuseStrength = 1;
props.SpecularColorReflectance = .5;
props.SpecularExponent = 20;
props.SpecularStrength = 1;
props.FaceColor= 'texture';
props.EdgeColor = 'none';
props.FaceLighting = 'phong';
props.Cdata = topo;

% Plot Earth
surface(x,y,z,props);
title('ECI Co-ordinates of Orbit')
axis equal
grid on
xlabel('x-axis (m)')
ylabel('y-axis (m)')
zlabel('z-axis (m)')
hold on

%plot ECI plot over 12H
for k = 1:31
plot3(R_ECI(1,:,k),R_ECI(2,:,k),R_ECI(3,:,k),'linewidth',1);
hold on
end
