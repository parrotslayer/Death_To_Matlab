%Testing Get_Star_Tracker
clc
clear all

deg2rad = pi/180;
r = 10e10;
Omega = 142.1145*deg2rad;
omega = 25.5213*deg2rad;
num_stars = 100;
inc = 98.7401*deg2rad;

%**************************** INPUT VARIABLES *************************
FOV = 90*deg2rad;   %field of view of star tracker
% Given ECEF position of the satellite and the Ground Station Coordinates
Sat_ECEF = [-3.9297e6;-4.5532e6;1.0266e7];
GS_LLH = [-33, 151,52];
%Convert satellite position to LLH
Sat_LLH = ECEF_to_LLH(Sat_ECEF(1),Sat_ECEF(2),Sat_ECEF(3));

t_epoch = [2003 03 18 0 0 0];
t_VE2002 = [2002 03 20 19 16 0];    %date of last VE UT
t_VE = etime(t_epoch,t_VE2002);


%get star positions (do outside loop?)
Star_ECI = Summon_Stars(r,inc,Omega,omega,num_stars);
visible_stars = Get_Star_Tracker(Sat_LLH,t_VE,Star_ECI,FOV);
