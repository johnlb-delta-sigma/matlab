%% dac_rz: Half-Return-to-Zero DAC's impulse response.
function [out] = dac_hz(t, t0,t1,t2, v1,v2)
	window = v1*(t>=t0 & t<=t1)-v2*(t>t1 & t<=t2);
	out = ( 1 )*window;
end