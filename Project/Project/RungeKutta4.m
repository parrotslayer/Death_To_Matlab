function [ xth ] = RungeKutta4(f, xt, ut, t, h, dynparams)
% Simple evaluation of 4th order Runge Kutta. Written for AMME5520 project
% by Ian Manchester.

    k1 = f(xt,ut,t, dynparams);
    k2 = f(xt+h/2*k1, ut, t+h/2, dynparams);
    k3 = f(xt+h/2*k2, ut, t+h/2, dynparams);
    k4 = f(xt+h*k3, ut, t+h, dynparams);
    
    xth = xt+h*(k1+2*k2+2*k3+k4)/6;

end

