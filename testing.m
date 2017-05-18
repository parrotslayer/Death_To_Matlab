%testing
clear all
close all 
clc

ephsat(1)= 7162*1000;               %a - semi major axis meters    
ephsat(2)=0.0000872;                %e - eccentricity deg
ephsat(3)=98.7401;                  %inc - inclination degrees
ephsat(4)=142.1145;                 %Omega - degrees
ephsat(5)=25.5213;                  %omega - degrees
ephsat(6)=283.570800000000;         %Mo - Mean Anomaly
ephsat(7)=2452716.5000;             %%Julian day = 2452716.5000 [days] 

[A,B,C] = Get_Estimated_ECEF(ephsat);

%% plots
figure
plot3(A(1,:),A(2,:),A(3,:),'*')
hold on
plot3(B(1,:),B(2,:),B(3,:),'.')
