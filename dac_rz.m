%% dac_rz: Return-to-Zero DAC's impulse response.
function [out] = dac_rz(t, t0,t1)
	window = (t>=t0 & t<=t1);
	out = ( 1 )*window;
end