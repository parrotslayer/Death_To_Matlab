function u = ComputeControl(u_0,currentstate,targetstate,modelparams)

K = LQRGain(modelparams);

deltaX = currentstate - targetstate;

u = u_0 - K * deltaX;