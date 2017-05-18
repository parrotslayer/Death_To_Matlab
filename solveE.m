function E = solveE(M,e,tol)
% Function solves Kepler's equation M = E-e*sin(E)
% Input - Mean anomaly M [rad] , Eccentricity e and Epsilon 
% Output  eccentric anomaly E [rad]. 
   	En0  = M;
    %use newtown's method to solve 
    %perform initial estimate
	En1 = En0 - (En0-e*sin(En0)- M)/(1 - e*cos(En0));
	%repeat until within tolerance
    while ( abs(En1-En0) > tol )
		En0 = En1;
		En1 = En0 - (En0 - e*sin(En0) - M)/(1 - e*cos(En0));
    end
	E = En1;

