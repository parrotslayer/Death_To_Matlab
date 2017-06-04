function Star_Constellation_ECI = Generate_Random_Stars(num_stars,radius)

rvals = 2*rand(num_stars,1)-1;
elevation = asin(rvals);
azimuth = 2*pi*rand(num_stars,1);
[Star_Constellation_ECI(:,1),Star_Constellation_ECI(:,2),Star_Constellation_ECI(:,3)] = sph2cart(azimuth,elevation,radius);

end