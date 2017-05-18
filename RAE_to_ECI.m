%% RAE to ECI
% uses functions to convert from RAE to ECI
% Takes in RAE input vector, the time and GS coordinates in LLH
% Uses RAE -> LGCV
% LGCV -> ECEF
% ECEF -> ECI
function output = RAE_to_ECI(input,time,GS_LLH)
LGCV = RAE_deg_to_LGCV(input);
ECEF = LGCV_to_ECEF(GS_LLH(1), GS_LLH(2), GS_LLH(3), LGCV(1), LGCV(2), LGCV(3));
output = ECEF_to_ECI([ECEF;time]);
end